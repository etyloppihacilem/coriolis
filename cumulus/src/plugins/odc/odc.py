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
from queue import LifoQueue
from threading import Event, Thread

from coriolis.Hurricane import Cell, Net
from sympy import S

from .CellODCCache import CellODCCache
from .FFDatabase import FFDatabase
from .ODCWalker import ODCWalker


class odc:
    def __init__(self, circuit: Cell, enable_estimate=True):
        self._cell: Cell = circuit
        self._todo: LifoQueue = LifoQueue()
        self._db: FFDatabase = FFDatabase(enable_estimate=enable_estimate)
        self._cache: CellODCCache = CellODCCache()
        self._done: Event = Event()
        self._enable_estimate = enable_estimate

        for net in self._cell.getExternalNets():
            if net.getDirection() == Net.Direction.OUT:
                self._todo.put(ODCWalker(net=net, todo=self._todo, results=self._db, cache=self._cache))

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

    def computeODC(self, force_simplify=False, refresh_rate=2, pretty=True):
        assert refresh_rate > 0
        print(f"Extracting observability for {self._cell.getName()}")
        started = datetime.now()
        print(f"Starting at {str(started).split('.')[0]}")
        runner = Thread(target=self.run_odc)
        runner.start()
        prev_walker_number = 0
        avg_walker_growth = []
        prev_walker_alive = 0
        avg_walker_alive = []
        prev_iter_count = 0
        avg_iter_speed = []
        print(f"Stats (live, refresh every {refresh_rate}s) :")
        while not self._done.is_set():
            print(f"Elapsed time : {str(datetime.now() - started).split('.')[0]}")
            walker_number = ODCWalker.walker_number
            walker_growth = (walker_number - prev_walker_number) / refresh_rate
            avg_walker_growth.append(walker_growth)
            print(f"Walker created : {walker_number} ({walker_growth:+} walker/s)")
            prev_walker_number = walker_number
            walker_alive = self._todo.qsize()
            avg_walker_alive.append(walker_alive)
            print(f"Walker alive : {walker_alive} ({(walker_alive - prev_walker_alive) / refresh_rate:+} walker/s)")
            prev_walker_alive = walker_alive
            iter_count = ODCWalker.iter_count
            iter_speed = (iter_count - prev_iter_count) / refresh_rate
            avg_iter_speed.append(iter_speed)
            print(f"Iterations : {iter_count} ({iter_speed:} iterations/s)")
            prev_iter_count = iter_count
            self._done.wait(timeout=refresh_rate)
            if pretty:
                print("\033[F\033[K" * 4, end="")
        if pretty:
            print("\033[F\033[K", end="")
        print("Simplifying functions, could take some time...")
        self._db.compute_functions(pretty)
        if pretty:
            print("\033[F\033[K", end="")
        runner.join()
        print("Stats :")
        print(f"  Elapsed time : {str(datetime.now() - started).split('.')[0]}")
        print(f"  Walkers : {ODCWalker.walker_number}")
        print(f"    avg. growth : {sum(avg_walker_growth) / max(len(avg_walker_growth), 1):+.2f} walker/s")
        print(f"    avg. alive  : {sum(avg_walker_alive) / max(len(avg_walker_alive), 1):.2f} walker")
        print(f"  Iterations : {ODCWalker.iter_count}")
        print(f"    avg. speed  : {sum(avg_iter_speed) / max(len(avg_iter_speed), 1):.2f} iteration/s")
        print(f"    it. per w.  : {ODCWalker.iter_count/max(ODCWalker.walker_number, 1):.2f} iteration/walker")
        print(f"  Results: {len(self._db)} flip-flops")
        functions = [f for f in self._db.values() if f.function != S.true]
        activation = len(functions)
        print(f"    With activation: {activation} flip-flops ({activation*100/max(len(self._db), 1):.2f}%)")
        print(f"    Simplified: {self._db.opti} functions out of {
              activation} ({self._db.opti*100/max(activation, 1):.2f}%)")
        print(f"    Variables removed: {self._db.variables_removed}")
        nb_var = [len(f.function.atoms()) for f in functions]
        print(f"    Avg. variables: {sum(nb_var)/max(len(nb_var), 1):.2f}")
        if len(nb_var) > 0:
            print(f"    Max. variables: {max(nb_var)}")
            print(f"    Min. variables: {min(nb_var)}")
        if self._enable_estimate:
            print(f"    Est. impact on cells: {self._db.compute_estimate()*100:.2f}%")
        # print("Iteration repartition")
        # for index, count in enumerate(ODCWalker.iter_rep):
        #     print(f"{index}: {count}")
        print("ODC done.")
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
        print(f"ODC results saved to {filename}")
