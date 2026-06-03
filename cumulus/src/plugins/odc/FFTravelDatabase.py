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
# |  Python      :   "./plugins/odc/FFTravelDatabase.py"            |
# +-----------------------------------------------------------------+

from coriolis.Hurricane import Cell

from .CellODC import CellODC


class FFTravelDatabase:
    def __init__(self, odc):
        self._ffs: set[str] = set()
        self.nets_true = set()
        self.odc = odc

    def __contains__(self, cell):
        if type(cell) is Cell:
            return cell.getName() in self._ffs
        return cell in self._ffs

    def addNewFF(self, ff: Cell, ff_info: CellODC, function, path, depth):
        ff_name = ff.getName()
        if ff_name in self._ffs:
            return True  # return true if walker should stop
        else:
            self._ffs.add(ff_name)
        return False

    def __len__(self):
        return len(self._ffs)
