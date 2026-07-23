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
# |  Python      :   "./plugins/odc/ODCGraph.py"                    |
# +-----------------------------------------------------------------+

from .InputsGraph import InputsGraph


class InputsGrapher:
    def __init__(self, info_cache):
        self.info_cache = info_cache
        self.node_db = {}
        self.input_nets = set()
        self.graphs = []

    def graph(self, instance):
        graph = InputsGraph(instance, self.info_cache, self.input_nets, self.node_db)
        self.graphs.append(graph)
