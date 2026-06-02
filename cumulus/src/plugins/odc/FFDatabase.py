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
# |  Python      :   "./plugins/odc/FFDatabase.py"                  |
# +-----------------------------------------------------------------+

from queue import Queue

from coriolis.Hurricane import Cell
from sympy import Not, Or, S, simplify_logic

from .CellODC import CellODC


def generateDepthName(ff, depth):
    if len(depth) > 0:
        return ff.getName() + "." + ".".join([c.getName() for c in depth])
    else:
        return ff.getName()


def optimize_FFEntry(entry):
    working = []
    already_opti = []
    while len(entry.functions):
        func = entry.functions.pop(0)
        entry.no_opti = Or(entry.no_opti, func)
        atoms = func.atoms() - {f.args[0] for f in func.atoms(Not)}
        atoms |= func.atoms(Not)
        working.append([atoms, func])
    working = sorted(working, key=lambda x: len(x[0]))
    while len(working) > 0:
        atoms, func = working.pop(0)
        done = False
        for i in range(len(working)):
            if atoms.issubset(working[i][0]):
                working[i][1] = simplify_logic(Or(func, working[i][1]), force=True)
                working[i][0] = working[i][1].atoms() - {
                    f.args[0] for f in working[i][1].atoms(Not)
                }
                working[i][0] |= working[i][1].atoms(Not)
                done = True
                break
        if not done:
            already_opti.append([atoms, func])
        else:
            working += already_opti
            already_opti.clear()
            working = sorted(working, key=lambda x: len(x[0]))
    working.clear()
    while len(already_opti) > 1:
        _, func = already_opti.pop(0)
        already_opti[0][1] = simplify_logic(Or(func, already_opti[0][1]))
    if len(already_opti) == 0:
        print(f"[ERROR] Simplifying function for {entry.name}, default to True.")
        print("")
        entry.function = S.true
        return
    entry.function = already_opti[0][1]


class FFEntry:
    def __init__(self):
        self.function: S = S.false
        self.no_opti: S = S.false
        self.functions = list()
        self.is_true = False
        self.name: str = ""

    def __str__(self):
        return f"{self.name}: {self.function}"


class FFDatabase:
    def __init__(self, odc):
        self._ff: dict[str, FFEntry] = {}
        self._ffs: set[str] = set()
        self._len = 0
        self.nets_true = set()
        self.opti = 0
        self.variables_removed = 0
        self.odc = odc
        self.path_ff: dict[int, set[str]] = {}

    def __contains__(self, cell):
        if type(cell) is Cell:
            return cell.getName() in self._ff
        return cell in self._ff

    def __getitem__(self, key):
        return self._ff[key]

    def addNewFF(self, ff: Cell, ff_info: CellODC, function, path, depth):
        ff_name = generateDepthName(ff, depth)
        if ff_name in self._ffs:
            self.path_ff[ff_name].update(path)
            old_entry = self._ff[ff_name]
            if old_entry.function == S.true or function == S.true:
                old_entry.function = S.true
                return True
            old_entry.functions.append(function)
            return True  # return true if walker should stop
        else:
            self._ffs.add(ff_name)
            entry = FFEntry()
            entry.function = function
            entry.functions.append(function)
            entry.name = ff_name
            self._ff[ff_name] = entry
            self.path_ff[ff_name] = set(path)
        return False

    def items(self):
        return self._ff.items()

    def values(self):
        return self._ff.values()

    def __len__(self):
        return len(self._ff)

    def compute_functions(self):
        from .odc import ODCVerbose

        for i, entry in enumerate(self._ff.values()):
            self.odc.printv(
                ODCVerbose.Normal,
                f"{i}/{len(self._ff)} ({i * 100 / len(self._ff):.2f}%), {
                    self.opti
                } optimizations found",
            )
            if entry.function == S.true:
                entry.functions.clear()
                entry.no_opti = S.true
                self.odc.erasev(ODCVerbose.Normal)
                continue
            optimize_FFEntry(entry)
            if entry.function != entry.no_opti:
                self.opti += 1
                self.variables_removed += len(entry.no_opti.atoms()) - len(
                    entry.function.atoms()
                )
            self.odc.erasev(ODCVerbose.Normal)

    def compute_estimate(self):
        results = {}
        for ff, cells in self.path_ff.items():
            affected = 1 if self._ff[ff].function != S.true else 0
            for cell in cells:
                try:
                    results[cell].append(affected)
                except KeyError:
                    results[cell] = [affected]
        calc_results = [sum(i) / len(i) for i in results.values()]
        return sum(calc_results) / len(calc_results)
