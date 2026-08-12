# This file is part of the Coriolis Software.
# Copyright (c) Sorbonne Université 2019-2026, All Rights Reserved
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

from coriolis.Hurricane import Cell
from sympy import S
from datetime import datetime

from .CellInfoCache import CellInfoCache
from .ODCGraph import ODCGraph
from .ODCGrapher import ODCGrapher


class ODCStats:
    def __init__(self):
        self.ff_count = 0
        self.functions = 0
        self.no_functions = 0
        self.false_functions = 0
        self.true_functions = 0
        self.children_count = []
        self.cut_count = []
        self.cut_points = []
        self.atoms_count = []
        self.cut_begin = None
        self.cut_end = None
        self.func_begin = None
        self.func_end = None

    def computeStats(self, grapher):
        for graph in grapher.graphs:
            self.add_graph(graph)

    def add_graph(self, graph):
        self.ff_count += 1
        if graph.function is None:
            self.no_functions += 1
        elif graph.function == S.false:
            self.false_functions += 1
        elif graph.function == S.true:
            self.true_functions += 1
        else:
            self.functions += 1
        cuts = graph.getCutSet()
        self.children_count.append(len(graph.top.children))
        self.cut_count.append(len(cuts))
        self.cut_points.extend([len(c) for c in cuts])
        self.atoms_count.append(
            len(graph.function.atoms()) if graph.function is not None else 0
        )

    def print(self):
        print("Stats:")
        print("------")
        print(f"FF count: {self.ff_count} ffs")
        if self.ff_count < 0:
            return
        mwidth = max(
            [
                len(str(x))
                for x in [
                    self.functions,
                    self.no_functions,
                    self.false_functions,
                    self.true_functions,
                ]
            ]
        )
        if self.ff_count > 0:
            print(
                f"With functions      : {self.functions:{mwidth}} ffs ({self.functions * 100 / self.ff_count:6.2f}% )"
            )
            value = self.no_functions
            percent = self.no_functions * 100 / self.ff_count
            print(f"With no functions   : {value:{mwidth}} ffs ({percent:6.2f}% )")
            value = self.false_functions
            percent = self.false_functions * 100 / self.ff_count
            print(f"With false functions: {value:{mwidth}} ffs ({percent:6.2f}% )")
            value = self.true_functions
            percent = self.true_functions * 100 / self.ff_count
            print(f"With true functions : {value:{mwidth}} ffs ({percent:6.2f}% )")
        if len(self.children_count) > 0:
            min_child = min(self.children_count)
            mean_child = sum(self.children_count) / len(self.children_count)
            max_child = max(self.children_count)
            print(f"Children count : {min_child} / {mean_child:.2f} / {max_child}")
        if len(self.cut_count) > 0:
            min_cut = min(self.cut_count)
            mean_cut = sum(self.cut_count) / len(self.cut_count)
            max_cut = max(self.cut_count)
            print(f"Cut count : {min_cut} / {mean_cut:.2f} / {max_cut}")
        if len(self.cut_points) > 0:
            min_cut = min(self.cut_points)
            mean_cut = sum(self.cut_points) / len(self.cut_points)
            max_cut = max(self.cut_points)
            print(f"Cut points : {min_cut} / {mean_cut:.2f} / {max_cut}")
        if len(self.atoms_count) > 0:
            min_atoms = min(self.atoms_count)
            mean_atoms = sum(self.atoms_count) / len(self.atoms_count)
            max_atoms = max(self.atoms_count)
            print(f"Atoms count : {min_atoms} / {mean_atoms:.2f} / {max_atoms}")
        # print(f"Cut time: {str(self.cut_end - self.cut_begin).split('.')[0]}")
        # print(f"Function time: {str(self.func_end - self.func_begin).split('.')[0]}")


class odc:
    def __init__(self, cell: Cell):
        self.cell = cell
        self.info_cache = CellInfoCache()
        self.grapher = ODCGrapher(self.info_cache)

    def computeODC(self):
        stats = ODCStats()
        print("Recherche des bascules")
        bascules = []
        for instance in self.cell.getInstances():
            cell_info = self.info_cache[instance]
            if cell_info.isFlipflop:
                bascules.append(instance)
        self.grapher.createGraphs(bascules)
        self.grapher.runAll(stats)
        print(self.grapher.spectralClustering())
        print("Fin, affichage des statistiques")
        stats.print()
