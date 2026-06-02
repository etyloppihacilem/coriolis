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
from enum import Enum
from queue import LifoQueue
from threading import Event, Thread

from coriolis.Hurricane import Cell, Net
from sympy import S

from .CellODCCache import CellODCCache
from .FFDatabase import FFDatabase
from .ODCWalker import ODCWalker, isHierarchical


def ODCselector(cell, top=True):
    first = True
    help = """Commandes:
  run nb [output_name]: run odc on circuit nb in list.
  list : list instances in current circuit.
  help : print this message.
  exit : exit the selector.
  back : go back one level
  open nb : open circuit nb."""
    if top:
        print(help)
    cells = [cell]
    for inst in cell.getInstances():
        if isHierarchical(inst):
            cells.append(inst.getMasterCell())
    while True:
        if first:
            for nb, c in enumerate(cells):
                if nb == 0:
                    print(f"[{nb}] {c.getName()} (current)")
                else:
                    print(f"[{nb}] {c.getName()}")
            first = False
        sel = input(": ")
        if sel == "exit":
            return True
        elif sel == "list":
            for nb, c in enumerate(cells):
                if nb == 0:
                    print(f"[{nb}] {c.getName()} (current)")
                else:
                    print(f"[{nb}] {c.getName()}")
        elif sel == "back":
            return False
        elif sel == "help":
            print(help)
        elif len(sel) >= len("open") and sel[0 : len("open")] == "open":
            if len(sel) <= len("open") + 1:
                print("Usage: run nb. See 'help' for details.")
                continue
            if int(sel[len("open") + 1 :]) == 0:
                print("Circuit 0 is the one currently opened.")
                continue
            try:
                if ODCselector(cells[int(sel[len("open") + 1 :])], False):
                    return True
            except IndexError:
                print(f"Circuit {sel} not in list. try 'list'.")
            for nb, c in enumerate(cells):
                if nb == 0:
                    print(f"[{nb}] {c.getName()} (current)")
                else:
                    print(f"[{nb}] {c.getName()}")
        elif len(sel) >= len("run") and sel[0 : len("run")] == "run":
            args = sel.split(" ")
            if len(args) < 2:
                print("Usage: run nb [output_name]. See 'help' for details.")
                continue
            try:
                odc_util = odc(cells[int(args[1])])
            except IndexError:
                print(f"Circuit {sel} not in list. try 'list'.")
                continue
            if len(args) < 3:
                output_name = f"{cells[int(args[1])].getName()}.odc"
            else:
                output_name = f"{args[2]}.odc"
            odc_util.computeODC()
            odc_util.save_to_file(filename=output_name)
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
        self._todo: LifoQueue = LifoQueue()
        self._db: FFDatabase = FFDatabase(self)
        self._cache: CellODCCache = CellODCCache()
        self._done: Event = Event()
        self.verbose = ODCVerbose(verbose, erase)

        for net in self._cell.getExternalNets():
            if net.getDirection() == Net.Direction.OUT:
                self._todo.put(
                    ODCWalker(
                        net=net, todo=self._todo, results=self._db, cache=self._cache
                    )
                )

    def printv(self, verbose, *args, **kwargs):
        if self.verbose.val >= verbose:
            print(*args, **kwargs)

    def erasev(self, verbose, lines=1):
        if self.verbose.erase() and self.verbose.val >= verbose:
            print("\033[F\033[K" * lines, end="")

    def run_odc(self):
        ODCWalker.walker_number = 0
        ODCWalker.iter_count = 0
        ODCWalker.iter_rep = list()
        try:
            while not self._todo.empty():
                task = self._todo.get()
                task.run()
        except BaseException as e:
            self._done.set()
            raise e
        self._done.set()

    def computeODC(self, force_simplify=False, refresh_rate=2):
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
        self.erasev(ODCVerbose.Normal, 2)
        self.printv(ODCVerbose.Normal, f"Started at {str(started).split('.')[0]}")
        self.printv(ODCVerbose.Mini, "Simplifying functions, could take some time...")
        self._db.compute_functions()
        self.erasev(ODCVerbose.Mini, 1)
        runner.join()
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

    def save_to_file(self, filename="odc_results.odc"):
        if not self._done.is_set():
            print("[ERROR] ODC was not calculated yet, can not save to file.")
        results = OrderedDict(sorted(self._db.items(), key=lambda item: item[0]))
        with open(filename, "w") as f:
            for value in results.values():
                f.write(f"{value.name}: {value.function}\n")
                # f.write(f"{value.name}: {value.no_opti}\n\n")
        self.printv(ODCVerbose.Normal, f"ODC results saved to {filename}")
