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

from sympy import And, S, Symbol, simplify_logic
from .HurricaneAsGraph import getParents, HurricaneHashasble
from coriolis.Hurricane import Net


class InputsNode:
    def __init__(self, graph, instance, is_top=False):
        self.instance = instance
        self.parents = []
        self.parents_net = {}
        self.children = []
        self.children_net = {}
        for parent, net in getParents(self.instance):
            if is_top and net not in graph.inputs:
                continue # nothing is saved as net is already defined before...
            if parent is None:
                graph.input_nets.add(HurricaneHashasble(net))
                continue
            parent_info = graph.info_cache[parent]
            if parent_info.isFlipflop:
                graph.input_nets.add(HurricaneHashasble(net))
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

    def computeFunction(self, pin_name: str = None):
        info = self.graph.info_cache[self.instance]
        symbol_map = getSymbolsMap(self.instance)
        function = None
        if pin_name is None:
            if len(info.pin_function) > 1:
                print(
                    f"[WARNING] Maybe using wrong function for cell {self.instance.getName()}"
                )
            function = list(info.pin_function.value())[0]
        else:
            if type(pin_name) is Net:
                for plug in self.instance.getPlugs():
                    if plug.getNet() == pin_name:
                        function = info.pin_function[pin_name.getName()]
                        break
            else:
                function = info.pin_function[pin_name]
        ext_expr = replaceSymbols(
            function,
            symbol_map,
        ) # fonction avec le nom des nets en symboles
        for parent, net in self.parents_net:
            pfunc = parent.computeFunction()
            ext_expr = ext_expr.subs(net, pfunc)
        return ext_expr
