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

from .ODCNode import ODCNode
from .Cuts import CutSetDB, Cut
from .InputsGrapher import InputsGrapher


class ODCGraph:
    def __init__(self, ff, info_cache):
        self.node_db = {}
        self.info_cache = info_cache
        self.top = ODCNode(self, ff, is_top=True)
        self.cell_info = info_cache[ff]
        if not self.cell_info.isFlipflop:
            print("[ERROR] Can not build graph from non-flipflop cell.")
            raise ValueError
        self.d_max = 20
        self.m = 11
        self.n_cap = 35
        self.max_inputs = 6
        self.cut_db = CutSetDB()

    def add_node(self, instance):
        try:
            return self.node_db[instance.getName()]
        except KeyError:
            node = ODCNode(self, instance)
            self.node_db[instance.getName()] = node
            return node

    def computeCuts(self):
        self.getCuts(self.top)
        # OPTI:
        # self.cut_db.clear_except(self.top)
        return self.getCutSet()

    def getCuts(self, node, level=0):
        if level > self.d_max or len(node.children) == 0:
            self.cut_db[node] = Cut(node)  # cut is added in a cut set automatically
            return
        for child in node.children:
            self.getCuts(child, level + 1)
        combined = self.cut_db[node.children[0]]
        for child in node.children[1:]:
            new_cuts = combined * self.cut_db[child]  # produit cartesien
            if new_cuts.maxSize() > self.m or len(new_cuts) > self.n_cap:
                self.cut_db[node] = Cut(node)  # added in cut set automatically
                return
            combined = new_cuts
        if level != 0:
            combined.add(Cut(node))
        self.cut_db[node] = combined

    def getCutSet(self):
        return self.cut_db[self.top]

    def computeFunctions(self):
        cuts = self.getCutSet()
        discarded = 0
        for cut in cuts:
            cut.computeCutGraph(self.top)
            cut.top.computeInputs()
            grapher = InputsGrapher(self.info_cache)
            for cnode in cut.top.input_nodes:
                grapher.graph(cnode.node.instance, cnode.inputs)
            print(len(grapher.input_nets))
            if len(grapher.input_nets) > self.max_inputs:
                discarded += 1
                continue
