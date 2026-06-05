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
# |  Python      :   "./plugins/odc/FastWalker.py"                  |
# +-----------------------------------------------------------------+

from queue import LifoQueue

from coriolis.Hurricane import Net, Plug

from .CellODCCache import CellODCCache
from .FFDatabase import generateDepthName
from .ODCWalker import isHierarchical


class FastWalker:
    def __init__(
        self,
        net: Net = None,
        todo: LifoQueue = None,
        cache: CellODCCache = None,
        other=None,
        from_plug=False,
        plug: Plug = None,
        net_db=None,
        cells_db=None,
        ffs=None,
    ):
        if net is None and plug is None:
            raise AttributeError
        if other is None:
            if cache is None or cells_db is None:
                print("[ERROR] Can not build walker without cache or results.")
                raise AttributeError
            self._cache = cache
            self._from_plug = from_plug
            self._net = net
            self._plug = plug
            self._todo = todo
            self._path = list()
            self._depth = list()
            self._net_db = net_db
            self._cells_db = cells_db
            self._ffs = ffs
        else:
            self._cache = other._cache  # not a copy
            self._from_plug = from_plug
            self._net = net
            self._plug = plug
            self._todo = other._todo  # not a copy
            self._path = list(other._path)
            self._depth = list(other._depth)
            self._net_db = other._net_db
            self._cells_db = other._cells_db
            self._ffs = other._ffs

    def fork(self, net: Net = None, plug: Plug = None):
        if net:
            self._todo.put(FastWalker(other=self, net=net))
        elif plug:
            self._todo.put(FastWalker(other=self, plug=plug, from_plug=True))
        else:
            print("[ERROR] Calling fork for ODCWalker with no args")
            raise AttributeError

    def iterate_over_net(self):
        self._plug = None
        if (
            len(self._depth) > 0
            and self._net.isExternal()
            and self._net.getDirection() == Net.Direction.IN
        ):
            # we 'teleport' out of hierarchical cell and continue.
            instance_frontiere = self._depth.pop()
            upper_plug = instance_frontiere.getPlug(self._net)
            if upper_plug:
                self._net = upper_plug.getNet()
            else:
                print("[WARNING] FastWalker could not get out of hierarchical cell.")
                return
        net = self._net
        first_plug = True
        for plug in net.getPlugs():
            master_net = plug.getMasterNet()
            if master_net.getDirection() != Net.Direction.OUT:
                continue
            if first_plug:
                self._plug = plug
                first_plug = False
            else:
                self.fork(plug=plug)

    def iterate_over_plug(self, instance, master_output):
        # Getting out of cell
        first = True
        self._net = None
        for plug in instance.getPlugs():  # should use list in odc_info
            net = plug.getNet()
            master_net = plug.getMasterNet()
            if master_net.isSupply() or master_net.isClock():
                continue
            if self._net_db is not None:
                self._net_db.append(
                    {
                        "depth": [i.getName() for i in self._depth],
                        "name": instance.getName(),
                        "io": master_net.getName(),
                        "net": net.getName(),
                    }
                )
            if master_net.getDirection() != Net.Direction.IN:
                continue

            if first:
                self._net = net
                first = False
            else:
                self.fork(net=net)

    def run(self):
        iter = 0
        instance = None
        while self._plug is not None or self._net is not None:
            iter += 1
            if not self._from_plug:
                self.iterate_over_net()
                if self._plug is None:  # no plugs were found
                    break
            else:
                self._from_plug = False
            # All plugs explored
            instance = self._plug.getInstance()
            if isHierarchical(instance):
                self._depth.append(instance)
                internalNet = self._plug.getMasterNet()
                self._net = internalNet
                self._plug = None
                continue  # we are traveling to internal net, not over plugs
            master_output = self._plug.getMasterNet()
            self._plug = None

            if instance.getName() in self._cells_db:
                break  # already visited this cell
            else:
                self._cells_db.add(instance.getName())
            if self._ffs is not None:
                odc_info = self._cache[instance]
                if odc_info.isFlipflop:
                    self._ffs.add(instance.getName())
            self._path.append(instance.getName())

            self.iterate_over_plug(instance, master_output)
