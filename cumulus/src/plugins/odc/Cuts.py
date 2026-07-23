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

from .CutNode import CutNode
from .ODCNode import ODCNode


class Cut:
    def __init__(self, construct=None):
        if construct is None:
            self.border = set()
            return
        if type(construct) is ODCNode:
            self.border = set()
            self.border.add(construct)
        else:
            self.border = set(construct)
        self.top = None
        self.max_inputs = 6

    def add(self, node: ODCNode):
        self.border.add(node)

    def __hash__(self):
        return hash(frozenset(self.border))

    def __or__(self, other):
        return Cut(self.border | other.border)

    def __len__(self):
        return len(self.border)

    def __iter__(self):
        return self.border.__iter__()

    def __repr__(self):
        return ", ".join([e.getName() for e in self])

    def __contains__(self, value):
        return value in self.border

    def computeCutGraph(self, top):
        self.top = CutNode(self, top)


class CutSet:
    def __init__(self, construct=None):
        if construct is None:
            self.cuts = set()
            return
        if type(construct) is ODCNode:
            self.cuts = set()
            self.cuts.add(Cut(construct))
        elif type(construct) is Cut:
            self.cuts = set()
            self.cuts.add(construct)
        else:
            self.cuts = set(construct)

    def add(self, cut: Cut):
        self.cuts.add(cut)

    def __mul__(self, other):
        new = set([c1 | c2 for c1 in self.cuts for c2 in other.cuts])
        return CutSet(new)

    def __or__(self, other):
        return CutSet(self.cuts | other.cuts)

    def maxSize(self):
        return max([len(c) for c in self.cuts])

    def __len__(self):
        return len(self.cuts)

    def __iter__(self):
        return self.cuts.__iter__()


class CutSetDB:
    def __init__(self):
        self.cut_sets = dict()

    def __getitem__(self, key: ODCNode) -> CutSet:
        try:
            return self.cut_sets[key]
        except KeyError:
            cut_set = CutSet()
            self.cut_sets[key] = cut_set
            return cut_set

    def __setitem__(self, key, value: CutSet):
        if type(value) is Cut:
            value = CutSet(value)
        self.cut_sets[key] = value

    def items(self):
        return self.cut_sets.items()

    def clear_except(self, top):
        tmp = self.cut_sets[top]
        self.cut_sets.clear()
        self.cut_sets[top] = tmp
