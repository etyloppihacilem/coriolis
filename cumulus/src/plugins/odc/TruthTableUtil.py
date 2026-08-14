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
# |  Python      :   "./plugins/odc/TruthTableUtil.py"              |
# +-----------------------------------------------------------------+

from dd.bdd import BDD

import sympy as sp


class BDDSymPyUtil:
    def __init__(self, expr, bdd_manager=None):
        """
        expr: Expression SymPy (ex: sp.And(a, sp.Or(b, c)))
        bdd_manager: Instance partagée de dd.BDD (facultatif)
        """
        self.expr = expr
        self.bdd = bdd_manager if bdd_manager is not None else BDD()

        self.symbols = sorted(list(self.expr.atoms(sp.Symbol)), key=lambda s: s.name)
        self.var_names = [s.name for s in self.symbols]

        for name in self.var_names:
            if name not in self.bdd.vars:
                self.bdd.declare(name)

        self.node = self._sympy_to_bdd(self.expr)

    def _sympy_to_bdd(self, expr):
        """Convertit récursivement un arbre d'expression SymPy en nœud BDD"""
        if expr is sp.true:
            return self.bdd.true
        if expr is sp.false:
            return self.bdd.false
        if isinstance(expr, sp.Symbol):
            return self.bdd.var(expr.name)

        if isinstance(expr, sp.Not):
            arg = self._sympy_to_bdd(expr.args[0])
            return self.bdd.apply('not', arg)

        if isinstance(expr, sp.And):
            res = self.bdd.true
            for arg in expr.args:
                res = res & self._sympy_to_bdd(arg)
            return res

        if isinstance(expr, sp.Or):
            res = self.bdd.false
            for arg in expr.args:
                res = res | self._sympy_to_bdd(arg)
            return res

        if isinstance(expr, sp.Xor):
            res = self._sympy_to_bdd(expr.args[0])
            for arg in expr.args[1:]:
                res = res ^ self._sympy_to_bdd(arg)
            return res

        raise TypeError(f"Type d'expression SymPy non pris en charge : {type(expr)}")

    def ONSetSize(self):
        if self.node == self.bdd.false or self.node == 0:
            return 0

        if self.node == self.bdd.true:
            return 1 << len(self.var_names) if self.var_names else 1

        if not self.var_names:
            return 0

        return self.bdd.count(self.node, nvars=len(self.var_names))

    def universal_quantification(self, variables_to_quantify):
        var_strs = {
            v.name if isinstance(v, sp.Symbol) else str(v)
            for v in variables_to_quantify
        }

        # Quantification universelle native en BDD
        # \forall_x f = f_{x=0} \cdot f_{x=1}
        quantified_node = self.node
        for var in var_strs:
            if var in self.bdd.vars:
                quantified_node = self.bdd.forall([var], quantified_node)

        # Reconstruire un objet avec le même BDD manager
        new_instance = BDDSymPyUtil.__new__(BDDSymPyUtil)
        new_instance.bdd = self.bdd
        new_instance.node = quantified_node
        new_instance.var_names = self.var_names
        new_instance.expr = None  # L'expression SymPy explicite devient optionnelle
        return new_instance

    def __and__(self, other):
        """ET logique entre deux BDDs."""
        new_instance = BDDSymPyUtil.__new__(BDDSymPyUtil)
        new_instance.bdd = self.bdd
        new_instance.node = self.node & other.node
        new_instance.var_names = sorted(
            list(set(self.var_names) | set(other.var_names))
        )
        new_instance.expr = None
        return new_instance

    def __or__(self, other):
        """OU logique entre deux BDDs."""
        new_instance = BDDSymPyUtil.__new__(BDDSymPyUtil)
        new_instance.bdd = self.bdd
        new_instance.node = self.node | other.node
        new_instance.var_names = sorted(
            list(set(self.var_names) | set(other.var_names))
        )
        new_instance.expr = None
        return new_instance
