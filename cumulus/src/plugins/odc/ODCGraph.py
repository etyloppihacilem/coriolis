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

from coriolis.Hurricane import Instance, Net

from .Cuts import Cut, CutSetDB
from .InputsGrapher import InputsGrapher
from .ODCNode import ODCNode


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


class ODCGraph:
    """
    Exactly one graph per flip-flop.
    """

    def __init__(self, ff: Instance, info_cache):
        self.node_db = {}
        self.info_cache = info_cache
        self.top = ODCNode(self, ff, is_top=True) # this is the flip flop.
        self.cell_info = info_cache[ff]
        if not self.cell_info.isFlipflop:
            print("[ERROR] Can not build graph from non-flipflop cell.")
            raise ValueError
        self.cut_db = CutSetDB()
        self.global_inputs: dict[ODCNode, dict[str, Net]] = InputDict() # inputs of the whole graph.

        # Parameters
        self.d_max = 20
        self.m = 11
        self.n_cap = 35
        self.max_inputs = 6

    def register_inputs(self, node):
        for plug in node.instance.getPlugs():
            master_net = plug.getMasterNet()
            if (
                master_net.getDirection() != Net.Direction.IN
                or master_net.isSupply()
                or master_net.isClock()
            ):
                continue
            net = plug.getNet()
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
        # self.cut_db.clear_except(self.top)
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


"""
Donc il faut une manière unifiée d'exprimer les graphes (parce que pendant la construction, on peut se retrouver à
parcourir à des étapes différentes la même node pour la parcourir pas du même sens et pas pour les mêmes raisons)
Donc dans l'ordre :
    - on fais le graphe depuis les cut et on enregistre toutes les inputs vides, puis on les strikes au fur et à mesure
    - on parcours ces inputs et on construit le graph dans l'autre sens en back track en enregistrant et strikant aussi
      les inputs que l'on crée sur le chemin
    --- là on vérifie si on a bien 6 entrées ou moins
    - on enregistre aussi les sorties primaires de tout le graphe normalement ce sont les cellules de cut...
    - Il y a aussi un problème dans la façon dont est construit le graph de cut: si on coupe une porte et qu'un de ses enfants est aussi dans la cut, le lien entre les deux est supprimé.
    - Donc il faut un graph global et rajouter dessus des vues qui sont des sortes de sous graphes qui se basent sur le graph de hurricane et l'agrandissent au besoin.
    CHANGEMENT D'APPROCHE

    On fait un seul graph qui représente tout hurricane et contient toutes les infos, notamment les fonctions exprimées par rapport aux nets, les enfants, les parents... Il faut un moyen de savoir si les parents et enfants ont étés explorés (certainement le tableau à none). On arrête d'utiliser une liste pour les parents et les enfants, maintenant c'est un dict qui fait [node, net], toutes les nodes sont du même type. on utilise que le nom des nets, tous les noms de pin devront être convertis. On affiche jamais la clock ni les alims. Le reset il faut vraiment y reflechir parce qu'on va pas clock gater le reset...
    Par défaut, l'exploration des nodes s'arrête à une bascule ou au bout du circuit. Elle peut être vers les enfants ou les parents. Elle se propage. Le graphe à un dict [nom_d'instance, node] (pour que ce soit hashable). On peut dire que si une node est déjà explorée dans un sens, on a pas besoin de le refaire. Mais attentions aux infos qui pourraient manquer si l'exploration n'est pas complète. Chaque node est du même objet. On fait des objets qui contiennent des données sur les nodes à la limite, on essaye de ne pas stocker de données dans les nodes parce que c'est super chiant a extraire après. par exemple pour la récursion dans le parcours du graphe...

    Pour les vues, c'est des surcharges des fonctions de parcours du graphe mais avec un paramètre de set qui contient des nodes. Pour tout parcours, on ajoute le set et si une node n'est pas dans le set, on considère qu'elle n'existe pas.

    Peut être arrêter de faire de la récursion et mettre en place des piles, ce sera plus propre.
"""
