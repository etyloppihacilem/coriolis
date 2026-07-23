# This file is part of the Coriolis Software.
# Copyright (c) Sorbonne Université 2019-2023, All Rights Reserved
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
        for pin, direction in self._pins_direction.items():
            print(f"#   {pin}:{direction}")
        print(f"# is FF        : {self._is_ff}")
