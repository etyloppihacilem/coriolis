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
# |  Python      :   "./plugins/odc/odc.py"                         |
# +-----------------------------------------------------------------+

from coriolis.Hurricane import Cell, Net


def getOutNets(instance: Cell):
    for plug in instance.getPlugs():
        net = plug.getNet()
        master_net = plug.getMasterNet()
        if (
            master_net.getDirection() != Net.Direction.OUT
            or master_net.isSupply()
            or master_net.isClock()
        ):
            continue
        yield net


def getInPlugs(net):
    for plug in net.getPlugs():
        master_net = plug.getMasterNet()
        if master_net.getDirection() != Net.Direction.IN:
            continue
        yield plug


def getChildren(instance: Cell):
    for net in getOutNets(instance):
        for plug in getInPlugs(net):
            yield plug.getInstance()
