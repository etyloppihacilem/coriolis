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

    def print(self, index=0, nl=False):
        print(index, end=" " if not nl else "\n")
        if index == 0:
            print(f"{index}: {self.function}")
        return [(index, self.function)]


class ODCFunctionGroup:
    def __init__(self):
        self.function = None
        self.bdd = None
        self.members = []

    def append(self, function: ODCFunction):
        import sympy as sp

        self.members.append(function)
        # TODO: Il faut absolument implémenter une meilleure simplification
        self.function = (
            sp.simplify(sp.And(self.function, function.function))
            if self.function is not None
            else function.function
        )
        self.bdd = BDDSymPyUtil(self.function)

    def print(self, index=0, nl=True):
        first = index == 0
        to_print = []
        print(index, end=" " if not nl else "\n")
        for index, member in enumerate(self.members):
            ret = member.print(index + 1, index == len(self.members) - 1)
            index += len(ret)
            to_print.extend(ret)
        if first:
            for index, function in to_print:
                print(f"{index}: {function}")
