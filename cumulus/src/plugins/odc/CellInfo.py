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
# |  Python      :   "./plugins/odc/CellInfo.py"                    |
# +-----------------------------------------------------------------+

import re

from coriolis.CRL import getLibertyGroupFromCell
from coriolis.Hurricane import Cell
from sympy import S, symbols, simplify_logic
from sympy.parsing.sympy_parser import (
    implicit_multiplication_application,
    parse_expr,
    standard_transformations,
)


def lib_expr_parsing(expression: str, pins: dict[str, S]):
    clean_expr = expression.replace("!", "~")
    clean_expr = clean_expr.replace('"', "")
    transformations = standard_transformations + (implicit_multiplication_application,)
    parsed_func = parse_expr(
        clean_expr, transformations=transformations, local_dict=pins
    )
    return simplify_logic(parsed_func)


class CellInfo:
    grp_extract = re.compile(r"\((.*)\)")

    def __init__(self, cell: Cell):
        try:
            master_cell = cell.getMasterCell()
        except AttributeError:
            master_cell = cell
        self._master_cell: str = master_cell.getName()
        self._pins_direction: dict[str, str] = {}
        self._is_ff: bool = False
        self._is_in_lib = True
        self.pin_function = {}

        grp = getLibertyGroupFromCell(master_cell)
        if grp is None:
            print(f"[WARNING] No library information attached to {self._master_cell}")
            self._is_in_lib = False
            return

        # adding pins to list
        pins = grp.getGroups("pin.*")
        for pin in pins:
            pin_name = CellInfo.grp_extract.search(pin.getGroupName()).group(1)
            pin_direction = pin.getAttribute("direction").getValue()
            self._pins_direction[pin_name] = pin_direction

        local_dict = {
            pin_name: symbols(pin_name) for pin_name in self._pins_direction.keys()
        }
        for pin in pins:
            pin_name = CellInfo.grp_extract.search(pin.getGroupName()).group(1)
            pin_direction = self._pins_direction[pin_name]
            if pin_direction == "output":
                function = pin.getAttribute("function")
                if function is None:
                    print(f"[WARNING] No function for output pin {pin_name} of {self._master_cell}")
                else:
                    function_str = function.getValue()
                    expr = lib_expr_parsing(function_str, local_dict)
                    self.pin_function[pin_name] = expr
        cg_integrated_cell = grp.getAttribute("clock_gating_integrated_cell")
        if cg_integrated_cell is not None:
            self._is_ff = True
            return

        # checking wether cell is a flip flop
        ff_grp = grp.getGroups("ff\\(.*")
        if (
            len(ff_grp) >= 1
        ):  # Attention, peut être ff et steering à la fois !!!! (reset...)
            ff_grp = ff_grp[0]
            # cell is ff
            self._is_ff = True
            # extracting outputs
            match = re.search(r"\((\w*),?(\w*).*\)", ff_grp.getGroupName())
            out = None
            inv_out = None
            if match:
                out, inv_out = match.groups()
            else:
                print("[ERROR] flip-flop groups is not correctly formatted.")
                raise SyntaxError
            # will only traver through inputs used in next_state parameter

    @property
    def isFlipflop(self):
        return self._is_ff

    def __bool__(self):
        return self._is_in_lib

    def print(self):
        print(f"# {self._master_cell}")
        print("# pins")
        try:
            for pin, direction in self._pins_direction.items():
                if direction == "output":
                    print(f"#   {pin}:{direction} ({self.pin_function[pin]})")
                else:
                    print(f"#   {pin}:{direction}")
        except KeyError:
            print(f"#   {self.pin_function}")
        print(f"# is FF        : {self._is_ff}")
