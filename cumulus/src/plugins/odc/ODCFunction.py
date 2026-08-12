# This file is part of the Coriolis Software.
# Copyright (c) Sorbonne Université 2019-2026, All Rights Reserved
#
# +-----------------------------------------------------------------+
# |                   C O R I O L I S                               |
# |      C u m u l u s  -  P y t h o n   T o o l s                  |
# |                     |
# |  Author      :                              Hippolyte MELICA    |
# |  E-mail      :   hippolyte.melica@etu.sorbonne-universite.fr    |
# | =============================================================== |
# |  Python      :   "./plugins/odc/ODCFunction.py"                 |
# +-----------------------------------------------------------------+

from .TruthTableUtil import BDDSymPyUtil


def similarity(r1, r2):
    f1 = r1.bdd
    f2 = r2.bdd

    vars1 = set(f1.var_names)
    vars2 = set(f2.var_names)
    common_vars = vars1 & vars2
    if len(common_vars) == 0:
        return 0.0

    z_vars = vars1 ^ vars2

    f_inter = f1 & f2
    size_inter = f_inter.ONSetSize()
    if size_inter == 0:
        return 0.0

    f_quantified = f_inter.universal_quantification(z_vars)
    size_quantified = f_quantified.ONSetSize()

    f_union = f1 | f2
    size_union = f_union.ONSetSize()
    if size_union == 0:
        return 0.0

    return (size_inter + size_quantified) / (2.0 * size_union)


class ODCFunction:
    def __init__(self, graph):
        self.info_cache = graph.info_cache
        self.top_net = graph.top_net
        self.instance = graph.top.link.instance  # because we don't use graph anymore
        self.function = graph.function
        self.bdd = BDDSymPyUtil(self.function)
