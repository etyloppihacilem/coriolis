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
from .HurricaneAsGraph import getChildren
from .Cuts import Cuts

d_max = 20
m = 11
n_cap = 35


def computeCuts(info_cache, result, ff, level=0):
    global d_max
    global m
    global n_cap
    if level > d_max:
        result[ff.getName()] = {Cuts({ff.getName()})}
        return
    children = [c for c in getChildren(ff)]
    if len(children) == 0:
        result[ff.getName()] = {Cuts({ff.getName()})}
        return
    for child in children:
        if info_cache[child].isFlipflop:
            continue # ERROR: il ne faut pas continuer, il faut marquer la sortie comme observée...
        computeCuts(info_cache, result, child, level + 1)
    combined = result[children[0].getName()]
    for child in children[1:]:
        new_cuts = set([Cuts(c1.union(c2)) for c1 in combined for c2 in result[child.getName()]])
        max_size = max([len(c) for c in new_cuts])
        if max_size > m or len(new_cuts) > n_cap:
            result[ff.getName()] = {Cuts({ff.getName()})}
            return
        combined = new_cuts
    result[ff.getName()] = combined.union({Cuts({ff.getName()})})


class odc:
    def __init__(self, cell: Cell):
        self.cell = cell
        self.info_cache = CellInfoCache()

    def computeODC(self):
        result = {}
        for instance in self.cell.getInstances():
            cell_info = self.info_cache[instance]
            if cell_info.isFlipflop:
                cuts = dict()
                computeCuts(self.info_cache, cuts, instance)
                result[instance.getName()] = cuts[instance.getName()]
        # import json
        # with open("odcgate_cut_test.json", "w") as f:
        #     json.dump(result, f, indent=2, default=lambda o: list(o))
        for inst, cuts in result.items():
            print(f"{inst} :")
            for cut in cuts:
                if type(cut) is Cuts:
                    print("  ---")
                    for c in cut:
                        print(f"    {c}")
                else:
                    print(f"  {type(cut)} -> {cut}")
