#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                                   Fahima Islam
#                           Oak Ridge National Laboratory
#
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


import numpy as np
from mcni.utils import conversion as conv
import unittest as unittest
from sample import Sample


class TestCase(unittest.TestCase):


    def test1(self):
        'Neutrons_from above'
        

        rq = lambda Q: np.exp(-Q*Q/25)

        s = Sample('sample', 5, 5, rq)

        from mcni import neutron_buffer, neutron
        N = 8
        nb = neutron_buffer(N)
        t = 1.
        for i in range(N):
            vi = np.array((0,-(i+1)*100,(i+1)*100))
            ri = -vi*t
            nb[i] = neutron(r=ri, v=vi, s=(0,0), time=0., prob=1.)
            continue

        nb2 = s.process(nb)
        for i, (n, n2) in enumerate(zip(nb, nb2)):
            # print "input:", n
            # print "output:", n2
            vf = np.array((0,(i+1)*100,(i+1)*100))
            np.allclose(vf, n2.state.velocity)
            np.allclose([0,0,0], n2.state.position)
            np.isclose(n2.time, 1.)
            vy = vf[1]
            Q = np.abs(vy*2) * conv.V2K
            np.isclose(rq(Q), n2.probability)
        return

    def test2(self):
        'Neutrons_from below'

        rq = lambda Q: Q * 0.5 + 1

        s = Sample('sample', 5, 5, rq)

        from mcni import neutron_buffer, neutron
        N = 8
        n = neutron_buffer(N)
        for i in range(N):
            n[i] = neutron(r=(0, -1, 0), v=(0, -3000, 0), s=(0, 1), time=1000., prob=10.)
            continue

        s.process(n)

        return

    def test3(self):
        'Neutrons parallel to the surface'

        rq = lambda Q: Q * 0.5 + 1

        s = Sample('sample', 5, 5, rq)

        from mcni import neutron_buffer, neutron
        N = 8
        n = neutron_buffer(N)
        for i in range(N):
            n[i] = neutron(r=(1, 0, 0), v=(3000, 3000, 0), s=(0, 1), time=1000., prob=10.)
            continue

        s.process(n)

        return

    def test4(self):
        'Neutrons passed the surface'

        rq = lambda Q: Q * 0.5 + 1

        s = Sample('sample', 5, 5, rq)

        from mcni import neutron_buffer, neutron
        N = 8
        n = neutron_buffer(N)
        for i in range(N):
            n[i] = neutron(r=(0, 0, 8), v=(0, 3000, 30000), s=(0, 1), time=1000., prob=10.)
            continue

        s.process(n)

        return
        

    pass # end of TestCase



if __name__ == "__main__": 
    unittest.main()

    
# End of file 
