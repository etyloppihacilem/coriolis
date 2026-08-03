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

import itertools

import numpy as np
from coriolis.Hurricane import Instance, Net
from sympy import lambdify, Or, And, simplify_logic, S

from .Cuts import Cut, CutSetDB, CutView
from .ODCNode import ODCNode


# The following is an implementation of the ODCGate algorithm presented in
# F. Liu, S. Yin, B. Jiang, Y. Ye, F. Farnia, and B. Yu,
# “ODCGate: Leveraging Observability Don’t Care Conditions for Enhanced Power Efficiency in Clock Gating,”
# ACM Trans. Des. Autom. Electron. Syst., p. 3812543, Apr. 2026, doi: 10.1145/3812543.


class InputDict:
    def __init__(self):
        self.data = {}

    def __setitem__(self, key, value):
        self.data[key] = value

    def __getitem__(self, key):
        try:
            return self.data[key]
        except KeyError:
            new_dict = dict()
            self.data[key] = new_dict
            return new_dict

    def __len__(self):
        return sum([len(value) for value in self.data.values()])

    def clean(self):
        keys_to_remove = [key for value, key in self.data.items() if len(value) == 0]
        for key in keys_to_remove:
            del self.data[key]

    def values(self):
        return self.data.values()

    def items(self):
        return self.data.items()

    def keys(self):
        return self.data.keys()


class ODCGraph:
    """
    Exactly one graph per flip-flop.
    """

    def __init__(self, ff: Instance, info_cache):
        self.excluded_net = {"rst_n"}  # TODO: should be used as parameter

        self.node_db = {}
        self.info_cache = info_cache
        self.top_net = None
        self.top = ODCNode(self, ff, is_top=True)  # this is the flip flop.
        self.cell_info = info_cache[ff]
        if not self.cell_info.isFlipflop:
            print("[ERROR] Can not build graph from non-flipflop cell.")
            raise ValueError
        self.cut_db = CutSetDB()
        # inputs of the whole graph.
        self.global_inputs: dict[ODCNode, dict[str, Net]] = InputDict()

        # Parameters
        self.d_max = 20
        self.m = 11
        self.n_cap = 35
        self.max_inputs = 6
        self.function = None

    def register_inputs(self, node):
        for plug in node.instance.getPlugs():
            master_net = plug.getMasterNet()
            if (
                # master_net.getDirection() != Net.Direction.IN
                master_net.isSupply() or master_net.isClock()
            ):
                continue
            net = plug.getNet()
            if master_net.getDirection() == Net.Direction.OUT:
                if node.is_top:
                    if self.top_net is not None:
                        print("[ERROR] Two output nets on top level cell.")
                    self.top_net = net
            elif master_net.getDirection() != Net.Direction.IN:
                continue
            if net.getName() in self.excluded_net:
                continue
            if not node.is_top:
                self.global_inputs[node][net.getName()] = net

    def connect(self, node, net):
        """
        Used when an instance wants to connect to a nodes input.
        """
        self.global_inputs[node].pop(net.getName(), None)

    def add_node(self, instance):
        try:
            return self.node_db[instance.getName()]
        except KeyError:
            node = ODCNode(self, instance)
            self.node_db[instance.getName()] = node
            return node

    def computeCuts(self):
        self.top.makeChildren()
        self.getCuts(self.top)
        # OPTI:
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
            print("[ERROR] Cut does not contain top_net.")
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
        to_expand = [node for node in self.global_inputs.keys()]
        for node in to_expand:
            node.makeParents()
        cuts = self.getCutSet()
        discarded = 0
        simulations = []
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

    def printStats(self):
        pass


"""
Donc il faut une manière unifiée d'exprimer les graphes (parce que pendant la construction, on peut se retrouver à
parcourir à des étapes différentes la même node pour la parcourir pas du même sens et pas pour les mêmes raisons)
Donc dans l'ordre :
  - on fais le graphe depuis les cut et on enregistre toutes les inputs vides, puis on les strikes au fur et à mesure
  - on parcours ces inputs et on construit le graph dans l'autre sens en back track en enregistrant et strikant aussi
    les inputs que l'on crée sur le chemin
  --- là on vérifie si on a bien 6 entrées ou moins
  - on enregistre aussi les sorties primaires de tout le graphe normalement ce sont les cellules de cut...
  - Il y a aussi un problème dans la façon dont est construit le graph de cut: si on coupe une porte et qu'un de ses
  enfants est aussi dans la cut, le lien entre les deux est supprimé.
  - Donc il faut un graph global et rajouter dessus des vues qui sont des sortes de sous graphes qui se basent sur le
  graph de hurricane et l'agrandissent au besoin.

CHANGEMENT D'APPROCHE

On fait un seul graph qui représente tout hurricane et contient toutes les infos, notamment les fonctions exprimées
par rapport aux nets, les enfants, les parents... Il faut un moyen de savoir si les parents et enfants ont étés
explorés (certainement le tableau à none). On arrête d'utiliser une liste pour les parents et les enfants,
maintenant c'est un dict qui fait [node, net], toutes les nodes sont du même type. on utilise que le nom des nets,
tous les noms de pin devront être convertis. On affiche jamais la clock ni les alims. Le reset il faut vraiment y
reflechir parce qu'on va pas clock gater le reset...
Par défaut, l'exploration des nodes s'arrête à une bascule ou au bout du circuit. Elle peut être vers les enfants
ou les parents. Elle se propage. Le graphe à un dict [nom_d'instance, node] (pour que ce soit hashable). On peut
dire que si une node est déjà explorée dans un sens, on a pas besoin de le refaire. Mais attentions aux infos qui
pourraient manquer si l'exploration n'est pas complète. Chaque node est du même objet. On fait des objets qui
contiennent des données sur les nodes à la limite, on essaye de ne pas stocker de données dans les nodes parce que
c'est super chiant a extraire après. par exemple pour la récursion dans le parcours du graphe...

Pour les vues, c'est des surcharges des fonctions de parcours du graphe mais avec un paramètre de set qui contient
des nodes. Pour tout parcours, on ajoute le set et si une node n'est pas dans le set, on considère qu'elle n'existe
pas.

Peut être arrêter de faire de la récursion et mettre en place des piles, ce sera plus propre.
"""
