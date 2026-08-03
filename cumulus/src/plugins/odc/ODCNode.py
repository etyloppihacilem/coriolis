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
from sympy import And, S, Symbol, simplify_logic


def getSymbolsMap(instance):
    ret = {}
    for plug in instance.getPlugs():
        net = plug.getNet()
        master_net = plug.getMasterNet()
        if net is None:
            continue
        ret[master_net.getName()] = Symbol(net.getName(), boolean=True)
    return ret


def replaceSymbols(expr, correspondance):
    ret = expr
    ret = ret.subs(correspondance)
    return ret


class ODCNode:
    def __init__(self, graph, instance: Instance, is_top=False):
        self.instance = instance
        self.graph = graph
        self.children = {}
        self.parents = {}
        # OPTI: is reset an input ?
        self.is_top = is_top
        self.graph.register_inputs(self)
        if not is_top:
            self.own_function = self.getOwnFunction()
            self.function = None
        else:
            self.function = Symbol(self.graph.top_net.getName())

    def getOwnFunction(self, pin_name: str = None):
        info = self.graph.info_cache[self.instance]
        symbol_map = getSymbolsMap(self.instance)
        function = None
        if pin_name is None:
            if len(info.pin_function) > 1:
                print(
                    f"[WARNING] More than one out pin. Maybe using wrong function for cell {self.instance.getName()}"
                )
            function = list(info.pin_function.values())[0]
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
        )  # fonction avec le nom des nets en symboles
        return ext_expr

    def computeFunction(self):
        if self.function is not None:
            return self.function
        ext_expr = self.own_function
        for parent, net in self.parents.items():
            pfunc = parent.computeFunction()
            ext_expr = ext_expr.subs(net.getName(), pfunc)
        self.function = ext_expr
        return ext_expr

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
            new_node.makeChildren()

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
            new_node.makeParents()

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

    def __repr__(self):
        return f"ODCNode({self.getName()}) {len(self.children)} children"
