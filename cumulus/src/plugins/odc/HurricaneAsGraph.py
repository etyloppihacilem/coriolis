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

from coriolis.Hurricane import Instance, Net


def getOutNets(instance: Instance):
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


def getInNets(instance: Instance):
    for plug in instance.getPlugs():
        net = plug.getNet()
        master_net = plug.getMasterNet()
        if (
            master_net.getDirection() != Net.Direction.IN
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


def getOutPlugs(net):
    for plug in net.getPlugs():
        master_net = plug.getMasterNet()
        if master_net.getDirection() != Net.Direction.OUT:
            continue
        yield plug


def getChildren(instance: Instance):
    for net in getOutNets(instance):
        found_plug = False
        for plug in getInPlugs(net):
            found_plug = True
            yield (plug.getInstance(), net)
        if not found_plug:
            yield (None, net)


cache = {}


def getParents(instance: Instance):
    global cache
    # amélioration : temps d'exécution /2
    try:
        return cache[instance.getName()]
    except KeyError:
        pass
    parents = []
    for net in getInNets(instance):
        found_plug = False
        for plug in getOutPlugs(net):
            found_plug = True
            parents.append((plug.getInstance(), net))
        if not found_plug:
            parents.append((None, net))
    cache[instance.getName()] = parents
    return parents


def printCache():
    global cache
    pass


class HurricaneHashasble:
    def __init__(self, obj):
        self.obj = obj

    def __hash__(self):
        return hash(self.obj.getName())
