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
# |  Python      :   "./plugins/odc/ODCGraph.py"                    |
# +-----------------------------------------------------------------+

import itertools

import numpy as np
from coriolis.Hurricane import Instance, Net
from sympy import lambdify, Or, And, simplify_logic, S

from .Cuts import Cut, CutSetDB, CutView
from .ODCNode import ODCNode
from .ODCGrapher import ODCGrapher


# The following is an implementation of the ODCGate algorithm presented in
# F. Liu, S. Yin, B. Jiang, Y. Ye, F. Farnia, and B. Yu,
# “ODCGate: Leveraging Observability Don’t Care Conditions for Enhanced Power Efficiency in Clock Gating,”
# ACM Trans. Des. Autom. Electron. Syst., p. 3812543, Apr. 2026, doi: 10.1145/3812543.


class ODCGraph:
    """
    Exactly one graph per flip-flop.
    """

    def __init__(self, ff: Instance, info_cache, grapher: ODCGrapher):
        self.excluded_net = {"rst_n"}  # TODO: should be used as parameter
        self.info_cache = info_cache

        if not self.info_cache[ff].isFlipflop:
            print("[ERROR] Can not build graph from non-flipflop cell.")
            raise ValueError
        self.grapher = grapher
        self.node_db = {}
        self.link_to_node = {}
        self.top_net = None
        self.top = ODCNode(self, ff, is_top=True)  # this is the flip flop.
        self.cut_db = CutSetDB()
        # inputs of the whole graph.
        self.graph_inputs: dict[ODCNode, Net] = {}
        self.function = None
        self.truth_table = None

        # Parameters
        # TODO: Should be set externally
        self.d_max = 20
        self.m = 11
        self.n_cap = 35
        self.max_inputs = 6

    def addNode(self, instance):
        try:
            return self.node_db[instance.getName()]
        except KeyError:
            pass
        node = ODCNode(self, instance)
        self.node_db[instance.getName()] = node
        return node

    def linkToNode(self, link):
        try:
            return self.link_to_node[link]
        except KeyError:
            pass
        return self.addNode(link.instance)

    def computeCuts(self):
        self.getCuts(self.top)
        # OPTI: décider si on garde ou pas pour l'économie de mémoire
        self.cut_db.clear_except(self.top)
        return self.getCutSet()

    def getCuts(self, node, level=0):
        if level > self.d_max or len(node.children) == 0:
            self.cut_db[node] = Cut(node)  # cut is added in a cut set automatically
            return
        for child in node.children:
            self.getCuts(child, level + 1)
        combined = self.cut_db[node.children_list[0]]
        for child in node.children_list[1:]:
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

    def getTruthTable(self, expr):
        all_symbols = sorted(list(expr.atoms()), key=lambda s: s.name)
        excluded_name = self.top_net.getName()
        target_symbols = [s for s in all_symbols if s.name != excluded_name]
        excluded_symbol = next(
            (s for s in all_symbols if s.name == excluded_name), None
        )
        if len(target_symbols) < 1:
            return (None, None)
        if excluded_symbol is None:
            # print(f"[ERROR] Cut for {self.top} has no top_net in expression.")
            return (None, None)
        num_vars = len(target_symbols)
        if num_vars > self.max_inputs:
            print(
                f"[WARNING] Cut does contain more than {self.max_inputs} and this should not be possible."
            )
        grid = np.array(list(itertools.product([0, 1], repeat=num_vars)), dtype=int)
        args = target_symbols + ([excluded_symbol] if excluded_symbol else [])
        func = lambdify(args, expr, modules="numpy")
        inputs = [grid[:, i] for i in range(num_vars)]
        inputs2 = [grid[:, i] for i in range(num_vars)]
        inputs.append(np.full(grid.shape[0], 0))
        inputs2.append(np.full(grid.shape[0], 1))
        raw_output = func(*inputs)
        raw_output2 = func(*inputs2)
        output = np.broadcast_to(raw_output, grid.shape[0]).astype(int)
        output2 = np.broadcast_to(raw_output2, grid.shape[0]).astype(int)
        result = (output == output2).astype(int)  # numpy XNOR because of dark int stuff
        truth_table = np.column_stack((grid, result))
        # # begin debug
        # import pandas as pd
        # header = [s.name for s in target_symbols] + ["XNOR"]
        # df = pd.DataFrame(truth_table, columns=header)
        # pd.set_option("display.max_rows", None)
        # print(df)
        # # end debug
        mask = truth_table[:, -1] == 1
        reduced_table = truth_table[mask]
        return (reduced_table, target_symbols)

    def computeFunctions(self):
        cuts = self.getCutSet()
        discarded = 0
        simulations = []
        no_top_net = 0
        for cut in cuts:
            view = CutView(self, cut)
            inputs = view.getInputs()
            if len(inputs) > self.max_inputs:
                discarded += 1
                continue
            # keeping this cut
            functions = view.getFunctions()
            for func in functions:
                sim = self.getTruthTable(func)
                if sim[0] is None:
                    no_top_net += 1
                    continue
                simulations.append(sim)
        # concatenating functions
        function = None
        for table, syms in simulations:
            minterms_expr = []
            for row in table:
                term_literal = []
                for val, sym in zip(row[:-1], syms):
                    if val == 1:
                        term_literal.append(sym)
                    else:
                        term_literal.append(~sym)
                minterms_expr.append(And(*term_literal))
            function = And(Or(*minterms_expr), function if function is not None else S.true)
        if function is None:
            return None
        self.function = simplify_logic(function)
        return self.function
