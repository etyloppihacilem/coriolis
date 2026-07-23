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
# |  Python      :   "./plugins/odc/CutNode.py"                     |
# +-----------------------------------------------------------------+

from .ODCNode import ODCNode
from .Cuts import Cut


class CutNode:
    def __init__(self, cut: Cut, node: ODCNode):
        self.node = node
        self.cut = cut
        self.end = self.node in self.cut
        self.inputs = self.node.inputs
        self.getName = self.node.getName
        self.children = [CutNode(self.cut, self.node.children)] if not self.end else []
        self.inputs = {}

    def computeInputs(self):
        for child in self.children:
            child.ComputeInputs()
            self.inputs |= child.inputs
        if not self.node.is_top:
            self.inputs |= self.node.inputs
