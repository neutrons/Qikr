import mcvine, mcvine.components
from mcni.AbstractComponent import AbstractComponent
from mcni.utils import conversion
import numpy as np
import os
from mcni import neutron_buffer, neutron

class Sample(AbstractComponent):

    def __init__(self, name, xwidth, zheight, rq):

        self.name = name
        self.xwidth = xwidth;
        self.zheight = zheight
        self.rq=rq
        return

    def process(self, neutrons): # input neutron is buffer
        if not len(neutrons):
            return
        from mcni.neutron_storage import neutrons_as_npyarr, ndblsperneutron # number of doubles per neutrons thats means each neutron is represented by x, y, z, vx, vy, vz, s1, s2, t, t0, p (10 double variables)

        arr = neutrons_as_npyarr(neutrons) #converting the input neutrons to array
        arr.shape = -1, ndblsperneutron
        x = arr[:, 0];
        y = arr[:, 1];
        z = arr[:, 2]
        vx = arr[:, 3];
        vy = arr[:, 4];
        vz = arr[:, 5]
        s1 = arr[:, 6];
        s2 = arr[:, 7];
        t = arr[:, 8];
        t0 = t.copy()
        p = arr[:, 9]

        # propagate to y = 0
        self._propagateToY0(x, y, z, vx, vy, vz, t)

        # Apply filter if area is positive
        assert self.xwidth > 0 and self.zheight > 0

        # Filter
        ftr = (x >= -self.xwidth / 2) * (x <= self.xwidth / 2) * (z >= -self.zheight / 2) * (z <= self.zheight / 2) * (
                    t > t0)

        # reflection
        vy[ftr] *= -1

        # adjust probability
        Q = np.abs(conversion.V2K * 2 * vy[ftr])
        p[ftr] *= self.reflectivity(self.rq, Q)

        Pneutrons = neutron_buffer(len(neutrons))
        Pneutrons.from_npyarr(arr)
        return  Pneutrons

    def _propagateToY0(self, x, y, z, vx, vy, vz, t):
        dt = -y / vy
        x += vx * dt
        y[:] = 0
        z+= vz*dt
        t += dt
        return

    def reflectivity(self, rq, *args, **kwargs):
        R = rq(*args, **kwargs)
        return R


