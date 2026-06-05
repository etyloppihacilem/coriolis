# This file is part of the Coriolis Software.
# Copyright (c) Sorbonne Université 2019-2023, All Rights Reserved
#
# +-----------------------------------------------------------------+
# |                   C O R I O L I S                               |
# |      C u m u l u s  -  P y t h o n   T o o l s                  |
# |                                                                 |
# |  Author      :                              Hippolyte MELICA    |
# |  E-mail      :   hippolyte.melica@etu.sorbonne-universite.fr    |
# | =============================================================== |
# |  Python      :   "./plugins/odc/odc.py"                         |
# +-----------------------------------------------------------------+

from collections import OrderedDict
from datetime import datetime
import json
from queue import LifoQueue
from threading import Event, Thread

from coriolis.Hurricane import Cell, Net
from sympy import S

from .CellODCCache import CellODCCache
from .FFDatabase import FFDatabase
from .FastWalker import FastWalker
from .ODCWalker import ODCWalker, isHierarchical


def ODCselector(cell, commands=[], top=None):
    if top is None:
        if len(commands) == 0:
            commands.insert(0, "list")
            commands.insert(0, "help")
        else:
            commands.insert(0, "list")
        top = odc(cell)
    else:
        commands.insert(0, "list")
    help = """Commandes:
  run [output_name] : run odc on top level circuit. Will write in file output_name_odc.json
  select nb : selects circuit nb for odc.
  reset : resets selection.
  list : list instances in current circuit.
  help : print this message.
  exit : exit the selector.
  back : go back one level
  nets : print nets.
  open nb : open circuit nb."""
    cells = []
    for inst in cell.getInstances():
        if isHierarchical(inst):
            cells.append(inst.getMasterCell())
    while True:
        if len(commands) > 0:
            sel = commands.pop(0)
            print(f": {sel}")
        else:
            sel = input(": ")
        if sel == "exit":
            return True
        elif sel == "list":
            for nb, c in enumerate(cells):
                print(f"[{nb}] {c.getName()}")
        elif sel == "back":
            return False
        elif sel == "help":
            print(help)
        elif sel == "reset":
            top.clear_selection()
        elif len(sel) >= len("nets") and sel[0 : len("nets")] == "nets":
            args = sel.split(" ")
            if len(args) < 2:
                output_name = f"{top.getName()}_correspondance.json"
            else:
                output_name = f"{args[1]}_correspondance.json"
            top.generate_net_to_plug(filename=output_name)
        elif len(sel) >= len("open") and sel[0 : len("open")] == "open":
            if len(sel) <= len("open") + 1:
                print("Usage: open nb. See 'help' for details.")
                continue
            try:
                if ODCselector(
                    cells[int(sel[len("open") + 1 :])], commands=commands, top=top
                ):
                    return True
            except IndexError:
                print(f"Circuit {sel} not in list. try 'list'.")
            except ValueError:
                print(f"'{sel[len('open') + 1]}' is not an index.")
            commands.insert(0, "list")
        elif len(sel) >= len("select") and sel[0 : len("select")] == "select":
            if len(sel) <= len("select") + 1:
                print("Usage: select nb. See 'help' for details.")
                continue
            try:
                top.select(cells[int(sel[len("select") + 1 :])])
            except IndexError:
                print(f"Circuit {sel} not in list. try 'list'.")
            except ValueError:
                print(f"'{sel[len('open') + 1]}' is not an index.")
        elif len(sel) >= len("run") and sel[0 : len("run")] == "run":
            args = sel.split(" ")
            if len(args) < 2:
                output_name = f"{top.getName()}_odc.json"
            else:
                output_name = f"{args[1]}_odc.json"
            top.computeODC()
            top.save_to_file(filename=output_name)
        else:
            print(f"Command '{sel}' not found. try 'help'.")


class ODCVerbose:
    No = 0
    Mini = 1
    Normal = 2
    Full = 3

    def __init__(self, val, erase: bool = True):
        self._erase = erase
        self.val = val

    def erase(self, new=None):
        if new is not None:
            self._erase = new
        return self._erase


class odc:
    def __init__(self, circuit: Cell, verbose=ODCVerbose.Normal, erase=True):
        self._cell: Cell = circuit
        self._selection: set[str] = set()
        self._todo: LifoQueue = LifoQueue()
        self._db: FFDatabase = FFDatabase(self, self._selection)
        self._cache: CellODCCache = CellODCCache()
        self._done: Event = Event()
        self._net_db = []
        self.nets_pairs = []
        self.verbose = ODCVerbose(verbose, erase)

    def printv(self, verbose, *args, **kwargs):
        if self.verbose.val >= verbose:
            print(*args, **kwargs)

    def erasev(self, verbose, lines=1):
        if self.verbose.erase() and self.verbose.val >= verbose:
            print("\033[F\033[K" * lines, end="")

    def clear_selection(self):
        self._selection.clear()

    def run_travel(self, cell, net_db, ffs_db, nets_pairs=None):
        self._done.clear()
        _cells_db = set()
        for net in cell.getExternalNets():
            if net.getDirection() == Net.Direction.OUT:
                self._todo.put(
                    FastWalker(
                        net=net,
                        todo=self._todo,
                        cells_db=_cells_db,
                        cache=self._cache,
                        net_db=net_db,
                        ffs=ffs_db,
                        nets_pairs=nets_pairs,
                    )
                )
        try:
            while not self._todo.empty():
                task = self._todo.get()
                task.run()
        except BaseException as e:
            self._done.set()
            raise e
        self._done.set()

    def select(self, cell: Cell, refresh_rate=2):
        assert refresh_rate > 0
        self.printv(ODCVerbose.Normal, f"Selecting flip-flops in {cell.getName()}")
        started = datetime.now()
        self.printv(ODCVerbose.Normal, f"Starting at {str(started).split('.')[0]}")
        selected_ffs = set()
        runner = Thread(target=self.run_travel, args=(cell, None, selected_ffs))
        runner.start()
        self.printv(ODCVerbose.Normal, f"Stats (live, refresh every {refresh_rate}s) :")
        while not self._done.is_set():
            self.printv(
                ODCVerbose.Mini,
                f"Elapsed time : {str(datetime.now() - started).split('.')[0]}",
            )
            self.printv(
                ODCVerbose.Normal,
                f"Selected ff : {len(selected_ffs)} ffs",
            )
            self._done.wait(timeout=refresh_rate)
            self.erasev(ODCVerbose.Mini, 1)
            self.erasev(ODCVerbose.Normal, 1)
        self._done.clear()
        self.erasev(ODCVerbose.Normal, 2)
        self.printv(ODCVerbose.Normal, f"Started at {str(started).split('.')[0]}")
        runner.join()
        self.printv(ODCVerbose.Mini, "Stats :")
        self.printv(
            ODCVerbose.Mini,
            f"  Elapsed time : {str(datetime.now() - started).split('.')[0]}",
        )
        self.printv(ODCVerbose.Mini, f"  Results: {len(selected_ffs)} flip-flops")
        self._selection |= selected_ffs
        self.printv(ODCVerbose.Mini, "Selection done.")

    def get_all_nets(self, refresh_rate=2):
        assert refresh_rate > 0
        self.printv(ODCVerbose.Normal, f"Correspondances of {self._cell.getName()}")
        started = datetime.now()
        self.printv(ODCVerbose.Normal, f"Starting at {str(started).split('.')[0]}")
        runner = Thread(
            target=self.run_travel,
            args=(self._cell, self._net_db, None, self.nets_pairs),
        )
        runner.start()
        self.printv(ODCVerbose.Normal, f"Stats (live, refresh every {refresh_rate}s) :")
        while not self._done.is_set():
            self.printv(
                ODCVerbose.Mini,
                f"Elapsed time : {str(datetime.now() - started).split('.')[0]}",
            )
            self.printv(
                ODCVerbose.Normal,
                f"Correspondances : {len(self._net_db)} nets",
            )
            self._done.wait(timeout=refresh_rate)
            self.erasev(ODCVerbose.Mini, 1)
            self.erasev(ODCVerbose.Normal, 1)
        self._done.clear()
        self.erasev(ODCVerbose.Normal, 2)
        self.printv(ODCVerbose.Normal, f"Started at {str(started).split('.')[0]}")
        runner.join()
        self.printv(ODCVerbose.Mini, "Stats :")
        self.printv(
            ODCVerbose.Mini,
            f"  Elapsed time : {str(datetime.now() - started).split('.')[0]}",
        )
        self.printv(ODCVerbose.Mini, f"  Results: {len(self._net_db)} correspondances")
        self.printv(ODCVerbose.Mini, f"  Nets pairs: {len(self.nets_pairs)} pairs")
        self.printv(ODCVerbose.Mini, "Correspondances done.")

    def run_odc(self):
        ODCWalker.walker_number = 0
        ODCWalker.iter_count = 0
        ODCWalker.iter_rep = list()
        self._done.clear()
        for net in self._cell.getExternalNets():
            if net.getDirection() == Net.Direction.OUT:
                self._todo.put(
                    ODCWalker(
                        net=net, todo=self._todo, results=self._db, cache=self._cache
                    )
                )
        try:
            while not self._todo.empty():
                task = self._todo.get()
                task.run()
        except BaseException as e:
            self._done.set()
            raise e
        self._done.set()

    def computeODC(self, refresh_rate=2):
        assert refresh_rate > 0
        self.printv(
            ODCVerbose.Normal, f"Extracting observability for {self._cell.getName()}"
        )
        started = datetime.now()
        self.printv(ODCVerbose.Normal, f"Starting at {str(started).split('.')[0]}")
        runner = Thread(target=self.run_odc)
        runner.start()
        prev_walker_number = 0
        avg_walker_growth = []
        prev_walker_alive = 0
        avg_walker_alive = []
        prev_iter_count = 0
        avg_iter_speed = []
        self.printv(ODCVerbose.Normal, f"Stats (live, refresh every {refresh_rate}s) :")
        while not self._done.is_set():
            self.printv(
                ODCVerbose.Mini,
                f"Elapsed time : {str(datetime.now() - started).split('.')[0]}",
            )
            walker_number = ODCWalker.walker_number
            walker_growth = (walker_number - prev_walker_number) / refresh_rate
            avg_walker_growth.append(walker_growth)
            self.printv(
                ODCVerbose.Full,
                f"Walker created : {walker_number} ({walker_growth:+} walker/s)",
            )
            prev_walker_number = walker_number
            walker_alive = self._todo.qsize()
            avg_walker_alive.append(walker_alive)
            self.printv(
                ODCVerbose.Full,
                f"Walker alive : {walker_alive} ({
                    (walker_alive - prev_walker_alive) / refresh_rate:+} walker/s)",
            )
            prev_walker_alive = walker_alive
            iter_count = ODCWalker.iter_count
            iter_speed = (iter_count - prev_iter_count) / refresh_rate
            avg_iter_speed.append(iter_speed)
            self.printv(
                ODCVerbose.Normal,
                f"Iterations : {iter_count} ({iter_speed:} iterations/s)",
            )
            prev_iter_count = iter_count
            self._done.wait(timeout=refresh_rate)
            self.erasev(ODCVerbose.Mini, 1)
            self.erasev(ODCVerbose.Full, 2)
            self.erasev(ODCVerbose.Normal, 1)
        self._done.clear()
        self.erasev(ODCVerbose.Normal, 2)
        self.printv(ODCVerbose.Normal, f"Started at {str(started).split('.')[0]}")
        runner.join()
        self.printv(ODCVerbose.Mini, "Simplifying functions, could take some time...")
        self._db.compute_functions()
        self.erasev(ODCVerbose.Mini, 1)
        self.printv(ODCVerbose.Mini, "Stats :")
        self.printv(
            ODCVerbose.Mini,
            f"  Elapsed time : {str(datetime.now() - started).split('.')[0]}",
        )
        self.printv(ODCVerbose.Full, f"  Walkers : {ODCWalker.walker_number}")
        self.printv(
            ODCVerbose.Full,
            f"    avg. growth : {
                sum(avg_walker_growth) / max(len(avg_walker_growth), 1):+.2f} walker/s",
        )
        self.printv(
            ODCVerbose.Full,
            f"    avg. alive  : {
                sum(avg_walker_alive) / max(len(avg_walker_alive), 1):.2f} walker",
        )
        self.printv(ODCVerbose.Normal, f"  Iterations : {ODCWalker.iter_count}")
        self.printv(
            ODCVerbose.Normal,
            f"    avg. speed  : {
                sum(avg_iter_speed) / max(len(avg_iter_speed), 1):.2f} iteration/s",
        )
        self.printv(
            ODCVerbose.Full,
            f"    it. per w.  : {
                ODCWalker.iter_count
                / max(ODCWalker.walker_number, 1):.2f} iteration/walker",
        )
        self.printv(ODCVerbose.Mini, f"  Results: {len(self._db)} flip-flops")
        selected = [f for f in self._db.values() if f.selected]
        selection = len(selected)
        self.printv(
            ODCVerbose.Mini,
            f"    Selected: {selection} flip-flops ({
                selection * 100 / max(len(self._db), 1):.2f}%)",
        )
        functions = [f for f in self._db.values() if f.function != S.true]
        activation = len(functions)
        self.printv(
            ODCVerbose.Mini,
            f"    With activation: {activation} flip-flops ({
                activation * 100 / max(len(self._db), 1):.2f}%)",
        )
        self.printv(
            ODCVerbose.Normal,
            f"    Simplified: {self._db.opti} functions out of {activation} ({
                self._db.opti * 100 / max(activation, 1):.2f}%)",
        )
        self.printv(
            ODCVerbose.Mini, f"    Variables removed: {self._db.variables_removed}"
        )
        nb_var = [len(f.function.atoms()) for f in functions]
        self.printv(
            ODCVerbose.Full,
            f"    Avg. variables: {sum(nb_var) / max(len(nb_var), 1):.2f}",
        )
        if len(nb_var) > 0:
            self.printv(ODCVerbose.Full, f"    Max. variables: {max(nb_var)}")
            self.printv(ODCVerbose.Full, f"    Min. variables: {min(nb_var)}")
        self.printv(
            ODCVerbose.Mini,
            f"    Est. impact on cells: {self._db.compute_estimate() * 100:.2f}%",
        )
        # self.printv(ODCVerbose.Mini, "Iteration repartition")
        # for index, count in enumerate(ODCWalker.iter_rep):
        #     self.printv(ODCVerbose.Mini, f"{index}: {count}")
        self.printv(ODCVerbose.Mini, "ODC done.")
        ODCWalker.walker_number = 0
        ODCWalker.iter_count = 0
        ODCWalker.iter_rep = list()

    def save_to_file(self, filename="odc_results.json"):
        results = OrderedDict(sorted(self._db.items(), key=lambda item: item[0]))
        with open(filename, "w") as f:
            f.write(
                json.dumps(
                    {
                        "circuit": self._cell.getName(),
                        "odc_results": [
                            v.as_dict()
                            for v in results.values()
                            if v.function != S.true
                        ],
                    }
                )
            )
        self.printv(ODCVerbose.Normal, f"ODC results saved to {filename}")

    def generate_net_to_plug(self, filename="nets_correspondance.json"):
        self.get_all_nets()
        nets = set()
        for entry in self._db.values():
            nets |= set([str(s) for s in entry.function.atoms()])
        group1 = set([i[0].getName() for i in self.nets_pairs])
        group2 = set([i[1].getName() for i in self.nets_pairs])
        group1_corr = {i[0].getName(): i[1].getName() for i in self.nets_pairs}
        group2_corr = {i[1].getName(): i[0].getName() for i in self.nets_pairs}
        to_add = set()
        for n in nets:
            if n in group1:
                to_add.add(group1_corr[n])
            elif n in group2:
                to_add.add(group2_corr[n])
        nets |= to_add
        selected = []
        for infos in self._net_db:
            if infos["net"] in nets:
                selected.append(infos)
        if len(selected) < len(nets):
            print("[WARNING] found less nets than there is variables.")
        with open(filename, "w") as f:
            f.write(json.dumps(selected, indent=2))
