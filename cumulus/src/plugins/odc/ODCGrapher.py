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
# |  Python      :   "./plugins/odc/ODCGrapher.py"                  |
# +-----------------------------------------------------------------+

from datetime import datetime

import numpy as np
from coriolis.Hurricane import Instance, Net
from sklearn.cluster import SpectralClustering
from sklearn.metrics import silhouette_score
from sympy import Symbol, S

from .CellInfoCache import CellInfoCache
from .HurricaneAsGraph import getChildren, getParents
from .ODCFunction import ODCFunction, ODCFunctionGroup, similarity


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


class NodeNet:
    def __init__(self, node, net):
        self.node = node
        self.link = node
        self.net = net

    def __hash__(self):
        return hash(self.net.getName() + self.node.getName())

    def __iter__(self):
        yield self.node
        yield self.net


class NodeLink:
    def __init__(self, grapher, instance: Instance):
        self.excluded_net = {"rst_n"}  # TODO: should be used as parameter
        self.grapher = grapher
        self.instance = instance
        self._children = {}
        self._children_done = False
        self._parents = {}
        self._parents_done = False
        self.own_function = self._getOwnFunction() if not self.isFF() else None
        # only contains unconnected inputs
        self.inputs: dict[str, Net] = {}
        self.primary_inputs = set()
        self.output = None
        self.isPrimaryOutput = False
        for plug in self.instance.getPlugs():
            master_net = plug.getMasterNet()
            if (
                # master_net.getDirection() != Net.Direction.IN
                master_net.isSupply() or master_net.isClock()
            ):
                continue
            net = plug.getNet()
            if master_net.getDirection() == Net.Direction.OUT:
                if self.output is not None:
                    print("[ERROR] Two output nets on top level cell.")
                self.output = net
                continue
            if master_net.getDirection() != Net.Direction.IN:
                continue
            if net.getName() in self.excluded_net:
                continue
            self.inputs[net.getName()] = net

    def getName(self):
        return self.instance.getName()

    def isFF(self):
        return self.grapher.info_cache[self.instance].isFlipflop

    @property
    def children(self):
        if not self._children_done:
            self._makeChildren()
            self._children_done = True
        return self._children

    @property
    def parents(self):
        if not self._parents_done:
            self._makeParents()
            self._parents_done = True
        return self._parents

    def _getOwnFunction(self, pin_name: str = None):
        info = self.grapher.info_cache[self.instance]
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

    def _addParent(self, parent, net):
        self._parents[parent] = net
        if not parent.isFF():  # connected to ff means primary input
            self.inputs.pop(net.getName(), None)

    def _addChild(self, child, net):
        if not child.isFF():
            self._children[child] = net
        if not self.isFF():  # connected to ff means primary input
            child.inputs.pop(net.getName(), None)

    def _makeChildren(self):
        """
        Please do not call that from public scope.
        """
        for child, net in getChildren(self.instance):
            if child is None:
                self.grapher.primary_outputs[net.getName()] = net
                self.isPrimaryOutput = True
                continue
            child_info = self.grapher.info_cache[child]
            if child_info.isFlipflop:
                self.grapher.primary_outputs[net.getName()] = net
                self.isPrimaryOutput = True
                continue
            new_link = self.grapher[child]
            new_link._addParent(self, net)
            self._children[new_link] = net
            # new_link.children  # generate value is does not exist...

    def _makeParents(self):
        for parent, net in getParents(self.instance):
            if parent is None:
                self.primary_inputs.add(NodeNet(self, net))
                continue
            parent_info = self.grapher.info_cache[parent]
            if parent_info.isFlipflop:
                self.primary_inputs.add(NodeNet(self, net))
                continue
            new_link = self.grapher[parent]
            new_link._addChild(self, net)
            self._parents[new_link] = net
            # new_link.parents  # generate value is does not exist...

    def getLeaves(self, ret=None):
        if ret is None:
            ret = set()
        if self in ret:
            return ret
        if self.isPrimaryOutput:
            ret.add(self)
        for child in self.children.keys():
            child.getLeaves(ret)
        return ret

    def getPrimaryInputs(self, traveled=None):
        ret = set()
        if traveled is None:
            traveled = set()
        if self in traveled:
            return ret
        traveled.add(self)
        ret |= self.primary_inputs
        if len(ret) > 6:
            return ret
        for parent in self.parents.keys():
            ret |= parent.getPrimaryInputs(traveled)
            if len(ret) > 6:
                return ret
        return ret

    def __hash__(self):
        return hash(self.instance.getName())


class ODCGrapher:
    def __init__(self, info_cache: CellInfoCache):
        self.database = {}
        self.info_cache = info_cache
        self.primary_outputs: dict[str, Net] = {}
        self.graphs = []
        self.results = []

    def clean(self):
        self.database.clear()
        self.primary_outputs.clear()
        self.graphs.clear()

    def clear(self):
        self.clean()
        self.results.clear()

    def __getitem__(self, key) -> NodeLink:
        try:
            return self.database[key.getName()]
        except KeyError:
            pass
        ret = NodeLink(self, key)
        self.database[key.getName()] = ret
        return ret

    def __iter__(self):
        return self.graphs.__iter__()

    def newGraph(self, instance):
        from .ODCGraph import ODCGraph

        graph = ODCGraph(instance, self.info_cache, self)
        self.graphs.append(graph)
        return graph

    def createGraphs(self, instances):
        print("Creation des graphes")
        from .ODCGraph import ODCGraph

        self.graphs = list([ODCGraph(i, self.info_cache, self) for i in instances])

    def runAll(self, stats=None):
        print("Calcul des coupes (prend un peu de temps)")
        cuts_begin = datetime.now()
        for graph in self.graphs:
            graph.computeCuts()
        cuts_end = datetime.now()
        print("Calcul des fonctions (peut prendre beaucoup de temps)")
        func_begin = datetime.now()
        for graph in self.graphs:
            graph.computeFunctions()
        func_end = datetime.now()
        stats.computeStats(self)
        for graph in self.graphs:
            if graph.function is None or graph.function == S.false:
                continue
            self.results.append(ODCFunction(graph))
        self.backup = list(self.results)
        self.clean()  # on supprime tout sauf les résultats
        print("Fusion des fonctions")
        # Affichage des temps
        print(f"Cuts done in {str(cuts_end - cuts_begin).split('.')[0]}")
        print(f"Func done in {str(func_end - func_begin).split('.')[0]}")
        total = (cuts_end - cuts_begin) + (func_end - func_begin)
        print(f"All done in {str(total).split('.')[0]}")

    def spectralClustering(self, k=8):  # TODO: trouver un moyen d'évaluer k
        # OPTI: Look at hMETIS tool for clustering
        N = len(self.results)
        affinity_matrix = np.zeros((N, N))

        for i in range(N):
            for j in range(N):
                if i == j:
                    affinity_matrix[i, j] = 1.0
                else:
                    affinity_matrix[i, j] = similarity(self.results[i], self.results[j])
        groups = {}
        scores = {}
        reverse_scores = {}
        for k in range(2, len(self.results) // 2):
            clustering = SpectralClustering(
                n_clusters=k,
                affinity="precomputed",
                assign_labels="kmeans",
                random_state=42,
            )
            labels = clustering.fit_predict(affinity_matrix)
            distance_matrix = 1.0 - affinity_matrix
            np.fill_diagonal(distance_matrix, 0.0)  # Sécurité pour éviter les -0.0

            # 2. Calcul du score de silhouette
            # On précise metric="precomputed"
            score = silhouette_score(distance_matrix, labels, metric="precomputed")
            scores[k] = score
            if score not in reverse_scores:
                reverse_scores[score] = k # on veut le k le plus petit
            groups[k] = labels
        max_score = max(list(scores.values()))
        min_score = min(list(scores.values()))
        mean_score = sum(list(scores.values())) / len(scores)

        min_k = reverse_scores[min_score]
        print(f"scores: {min_score} / {mean_score} / {max_score}")
        print(f"Optimal k value is {min_k}")
        print(f"{len(scores)} vs. {len(reverse_scores)}")
        print("k,score")
        for k, score in scores.items():
            print(f"{k},{score}")
        return (groups[min_k], min_k)

    def runSpectralClustering(self):
        # TODO: Trouver une condition d'arrêt
        # En gros le débat pour la condition d'arrêt c'est qu'il semblerait qu'on ne fasse
        # du clustering que sur un niveau dans le papier. Ce serait intéressant de regarder
        # comment ça se comporte lorsqu'on fait un arbre de clusters.
        while len(self.results) > 1:
            print(len(self.results))
            groups, score = self.spectralClustering()
            new_results = []
            group_indexes = {}
            for index, g in enumerate(groups):
                group_index = None
                try:
                    group_index = group_indexes[g]
                except KeyError:
                    group_index = len(new_results)
                    group_indexes[g] = group_index
                    new_results.append(ODCFunctionGroup())
                new_results[group_index].append(self.results[index])
            self.results = new_results
            break
            # la recherche de k optimal rend l'itération ici inutile.

    def reset_results(self):
        self.results = list(self.backup)
