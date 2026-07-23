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
from .HurricaneAsGraph import getChildren


class ODCNode:
    def __init__(self, graph, instance: Instance, is_top=False):
        self.instance = instance
        self.info = graph.info_cache[instance]
        self.children = []
        self.children_net = {}
        self.parents = []
        self.parents_net = {}
        for child, net in getChildren(self.instance):
            if child is None:
                continue
            child_info = graph.info_cache[child]
            if child_info.isFlipflop:
                continue
            new_node = graph.add_node(child)
            new_node.markInput(self, net)
            self.children.append(new_node)
            self.children_net[new_node] = net
        self.inputs: dict[str, Net] = {}  # clock or power supply are no input
        # OPTI: is reset an input ?
        self.is_top = is_top
        if not is_top:
            for plug in self.instance.getPlugs():
                master_net = plug.getMasterNet()
                if (
                    master_net.getDirection() != Net.Direction.IN
                    or master_net.isSupply()
                    or master_net.isClock()
                ):
                    continue
                net = plug.getNet()
                self.inputs[net.getName()] = net

    def markInput(self, parent, net):
        self.inputs.pop(net.getName(), None)
        self.parents.append(parent)
        self.parents_net[parent] = net

    def __hash__(self):
        return hash(self.instance.getName())

    def getName(self):
        return self.instance.getName()
