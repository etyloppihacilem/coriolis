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
            print(
                f"With no functions   : {self.no_functions:{mwidth}} ffs ({
                    self.no_functions * 100 / self.ff_count:6.2f}% )"
            )
            print(
                f"With false functions: {self.false_functions:{mwidth}} ffs ({
                    self.false_functions * 100 / self.ff_count:6.2f}% )"
            )
            print(
                f"With true functions : {self.true_functions:{mwidth}} ffs ({
                    self.true_functions * 100 / self.ff_count:6.2f}% )"
            )
        if len(self.children_count) > 0:
            print(
                f"Children count : {min(self.children_count)} / {
                    sum(self.children_count) / len(self.children_count):.2f} / {
                    max(self.children_count)
                }"
            )
        if len(self.cut_count) > 0:
            print(
                f"Cut count : {min(self.cut_count)} / {
                    sum(self.cut_count) / len(self.cut_count):.2f} / {
                    max(self.cut_count)
                }"
            )
        if len(self.cut_points) > 0:
            print(
                f"Cut points : {min(self.cut_points)} / {
                    sum(self.cut_points) / len(self.cut_points):.2f} / {
                    max(self.cut_points)
                }"
            )
        if len(self.atoms_count) > 0:
            print(
                f"Atoms count : {min(self.atoms_count)} / {
                    sum(self.atoms_count) / len(self.atoms_count):.2f} / {
                    max(self.atoms_count)
                }"
            )
        print(f"Cut time: {str(self.cut_end - self.cut_begin).split('.')[0]}")
        print(f"Function time: {str(self.func_end - self.func_begin).split('.')[0]}")
        print(f"Graph creation time: {str(ODCGraph.creation_time).split('.')[0]}")
        print(f"Graph extension time: {str(ODCGraph.total_time).split('.')[0]}")


class odc:
    def __init__(self, cell: Cell):
        self.cell = cell
        self.info_cache = CellInfoCache()
        self.grapher = ODCGrapher(self.info_cache)

    def computeODC(self):
        stats = ODCStats()
        result = []
        print("Recherche des bascules")
        bascules = []
        for instance in self.cell.getInstances():
            cell_info = self.info_cache[instance]
            if cell_info.isFlipflop:
                bascules.append(instance)
        print("Début du calcul des cuts")
        stats.cut_begin = datetime.now()
        for index, instance in enumerate(bascules):
            print(f"Calculating cut {index} / {len(bascules)}")
            graph = ODCGraph(instance, self.info_cache, self.grapher)
            graph.computeCuts()
            result.append(graph)
        stats.cut_end = datetime.now()
        # import json
        # with open("odcgate_cut_test.json", "w") as f:
        #     json.dump(result, f, indent=2, default=lambda o: list(o))
        print("Calcul des fonctions")
        stats.func_begin = datetime.now()
        for index, graph in enumerate(result):
            # print(f"\n{graph.top.getName()} : {graph.top_net.getName()}")
            # function = graph.computeFunctions()
            print(f"Calcul de {index} / {len(result)}")
            graph.computeFunctions()
            # print(function)
            stats.add_graph(graph)
        stats.func_end = datetime.now()
        print("Fin, affichage des statistiques")
        stats.print()
