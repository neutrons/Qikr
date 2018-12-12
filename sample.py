import numpy as np
import mcvine, mcvine.components
from mcni.AbstractComponent import AbstractComponent
from mcni.utils import conversion
from mcni import neutron_buffer, neutron
from mcni.neutron_storage import neutrons_as_npyarr, ndblsperneutron

DEBUG = False


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

        # Number of doubles per neutrons thats means each neutron is represented
        # by x, y, z, vx, vy, vz, s1, s2, t, t0, p (10 double variables)
        arr = neutrons_as_npyarr(neutrons)
        arr.shape = -1, ndblsperneutron
        x = arr[:, 0];
        y = arr[:, 1];
        z = arr[:, 2]
        vx = arr[:, 3];
        vy = arr[:, 4];
        vz = arr[:, 5]
        t = arr[:, 8];
        t0 = t.copy()
        p = arr[:, 9]

        # Propagate to y = 0
        self._propagateToY0(x, y, z, vx, vy, vz, t)

        # Apply filter if area is positive
        assert self.xwidth > 0 and self.zheight > 0

        # Filter
        do_filter = True
        ftr = (x >= -self.xwidth / 2) * (x <= self.xwidth / 2) * (y >= -self.zheight / 2) * (y <= self.zheight / 2) * (t > t0)

        # Reflection
        speed = np.sqrt(vx**2+vy**2+vz**2)
        if do_filter:
            vy[ftr] *= -1
            #q = 2.0*conversion.V2K*speed[ftr] * vy[ftr]/vz[ftr]
            q = 2.0*conversion.V2K*vy[ftr]
            p[ftr] *= self.reflectivity(self.rq, q)
        else:
            vy *= -1
            #q = 2.0*conversion.V2K*speed * vy/vz
            q = 2.0*conversion.V2K*vy
            p *= self.reflectivity(self.rq, q)

        if DEBUG:
            for i in range(len(arr)):
                self.print_neutron(arr[i], q[i])

        neutrons.from_npyarr(arr)
        return  neutrons

    def print_neutron(self, n, q0=0.0):
        tof = n[8]
        vx = n[3]
        vy = n[4]
        vz = n[5]
        h = 6.626e-34  # m^2 kg s^-1
        m = 1.675e-27  # kg
        dist = tof * vz
        cst = dist / h * m
        wl  = tof / cst * 1e10
        
        speed = np.sqrt(vx**2+vy**2+vz**2)
        k = conversion.V2K*speed
        wl2 = 2.0 * np.pi / k


        q = 4.0*np.pi/wl * vy/vz
        theta = np.arcsin(vy/speed) * 180.0 / np.pi
        print "z = %s    wl = %s (%s)   q = %s (%s)   theta = %s" % (dist, wl, wl2, q, q0, theta)

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
