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
# |  Python      :   "./plugins/odc/ODCNode.py"                     |
# +-----------------------------------------------------------------+

from coriolis.Hurricane import Instance, Net
from .HurricaneAsGraph import getChildren, getParents, HurricaneHashasble


class ODCNode:
    def __init__(self, graph, instance: Instance, is_top=False):
        self.instance = instance
        self.graph = graph
        self.children = {}
        self.parents = {}
        # OPTI: is reset an input ?
        self.is_top = is_top
        if not is_top:
            self.graph.register_inputs(self)

    def makeChildren(self):
        for child, net in getChildren(self.instance):
            if child is None:
                continue
            child_info = self.graph.info_cache[child]
            if child_info.isFlipflop:
                continue
            new_node = self.graph.add_node(child)
            new_node.addParent(self, net)
            self.graph.connect(new_node, net)
            self.children[new_node] = net

    def makeParents(self):
        for parent, net in getParents(self.instance):
            if parent is None:
                continue
            parent_info = self.graph.info_cache[parent]
            if parent_info.isFlipflop:
                continue
            new_node = self.graph.add_node(parent)
            new_node.addChild(self, net)
            self.graph.connect(self, net)
            self.parents[new_node] = net

    def addParent(self, parent, net):
        self.parents[parent] = net

    def addChild(self, child, net):
        self.children[child] = net

    def __hash__(self):
        return hash(self.instance.getName())

    def getName(self):
        return self.instance.getName()

    @property
    def children_list(self):
        return list(self.children.keys())

    @property
    def parent_list(self):
        return list(self.children.keys())
