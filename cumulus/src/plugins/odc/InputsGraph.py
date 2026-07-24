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
# |  Python      :   "./plugins/odc/InputsGraph.py"                 |
# +-----------------------------------------------------------------+

from .InputsNode import InputsNode


class InputsGraph:
    def __init__(self, instance, inputs, info_cache, input_nets, node_db):
        self.info_cache = info_cache
        self.input_nets = input_nets
        self.node_db = node_db
        self.inputs = inputs # this contains nets that are primary inputs of the cut
        self.top = self.add_node(instance, is_top=True)

    def add_node(self, instance, is_top=False):
        try:
            return self.node_db[instance.getName()]
        except KeyError:
            node = InputsNode(self, instance, is_top)
            self.node_db[instance.getName()] = node
            return node

    def computeFunction(self):
        return self.top.computeFunction()
