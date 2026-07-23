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
    def __init__(self, instance, info_cache):
        self.info_cache = info_cache
        self.inputs_nets = set()
        self.node_db = set()
        self.top = self.add_node(instance)

    def add_node(self, instance):
        try:
            return self.node_db[instance.getName()]
        except KeyError:
            node = InputsNode(self, instance)
            self.node_db[instance.getName()] = node
            return node
