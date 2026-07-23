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
# |  Python      :   "./plugins/odc/InputsNode.py"                  |
# +-----------------------------------------------------------------+

from .HurricaneAsGraph import getParents, HurricaneHashasble


class InputsNode:
    def __init__(self, graph, instance):
        self.instance = instance
        self.parents = []
        self.parents_net = {}
        self.children = []
        self.children_net = {}
        for parent, net in getParents(self.instance):
            if parent is None:
                graph.inputs_nets.add(HurricaneHashasble(net))
                continue
            parent_info = graph.info_cache[parent]
            if parent_info.isFlipflop:
                graph.inputs_nets.add(HurricaneHashasble(net))
                continue
            new_node = graph.add_node(parent)
            self.parents.append(new_node)
            self.parents_net[new_node] = net
            new_node.markOutput(self, net)

    def __hash__(self):
        return hash(self.instance.getName())

    def markOutput(self, child, net):
        self.children.append(child)
        self.children_net[child] = net
