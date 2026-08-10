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
from .ODCGrapher import NodeNet
from sympy import And, S, Symbol, simplify_logic, Function


class ODCNode:
    def __init__(self, graph, instance: Instance, is_top=False):
        self.is_top = is_top
        self.link = graph.grapher[instance]
        self.graph = graph
        self.graph.link_to_node[self.link] = self
        self.inputs = None
        if is_top:
            self.graph.top_net = self.link.output
            self.function = Function("")(Symbol(self.graph.top_net.getName()))
        else:
            self.function = None

    @property
    def children(self):
        return {self.graph.linkToNode(link): net for link, net in self.link.children.items() if not link.isFF()}

    @property
    def parents(self):
        return {self.graph.linkToNode(link): net for link, net in self.link.parents.items() if not link.isFF()}

    def computeFunction(self):
        if self.function is not None:
            return self.function # will handle case where is_top is True
        ext_expr = self.link.own_function
        for parent, net in self.parents.items():
            pfunc = parent.computeFunction()
            ext_expr = ext_expr.subs(net.getName(), pfunc)
        self.function = ext_expr
        return ext_expr

    def getGraphInputs(self):
        if self.inputs is not None:
            return self.inputs
        self.inputs: set[NodeNet] = set()
        leaves = self.link.getLeaves()
        for leaf in leaves:
            for node, net in leaf.getPrimaryInputs():
                if net != self.graph.top_net:
                    self.inputs.add((net, node))
                if len(self.inputs) > self.graph.max_inputs:
                    return self.inputs
        return self.inputs

    def __hash__(self):
        return self.link.__hash__()

    def getName(self):
        return self.link.getName()

    @property
    def children_list(self):
        return list(self.children.keys())

    @property
    def parent_list(self):
        return list(self.children.keys())

    def __repr__(self):
        return f"ODCNode({self.getName()}) {len(self.children)} children"
