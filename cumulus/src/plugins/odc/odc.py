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

from .CellInfoCache import CellInfoCache
from .ODCGraph import ODCGraph


class odc:
    def __init__(self, cell: Cell):
        self.cell = cell
        self.info_cache = CellInfoCache()

    def computeODC(self):
        result = []
        for instance in self.cell.getInstances():
            cell_info = self.info_cache[instance]
            if cell_info.isFlipflop:
                graph = ODCGraph(instance, self.info_cache)
                graph.computeCuts()
                result.append(graph)
        # import json
        # with open("odcgate_cut_test.json", "w") as f:
        #     json.dump(result, f, indent=2, default=lambda o: list(o))
        for graph in result:
            print(f"{graph.top.getName()} :")
            for cut in graph.getCutSet():
                print(f"  {cut}")
