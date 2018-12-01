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


interactive = False


import unittest as unittest


class TestCase(unittest.TestCase):


    def test1(self):
        'Neutrons_from above'
        
        from sample_component import Sample

        rq = lambda Q: Q * 0.5 + 1

        s = Sample('sample', 5, 5, rq)

        from mcni import neutron_buffer, neutron
        N = 8
        n = neutron_buffer(N)
        for i in range(N):
            n[i] = neutron(r=(0,1,0), v=(0,3000,0), s=(0,1), time=1000., prob=10.)
            continue

        s.process(n)

        r = s.events

        return

    def test2(self):
        'Neutrons_from below'

        from sample_component import Sample

        rq = lambda Q: Q * 0.5 + 1

        s = Sample('sample', 5, 5, rq)

        from mcni import neutron_buffer, neutron
        N = 8
        n = neutron_buffer(N)
        for i in range(N):
            n[i] = neutron(r=(0, -1, 0), v=(0, -3000, 0), s=(0, 1), time=1000., prob=10.)
            continue

        s.process(n)

        r = s.events
        return

    def test3(self):
        'Neutrons parallel to the surface'

        from sample_component import Sample

        rq = lambda Q: Q * 0.5 + 1

        s = Sample('sample', 5, 5, rq)

        from mcni import neutron_buffer, neutron
        N = 8
        n = neutron_buffer(N)
        for i in range(N):
            n[i] = neutron(r=(1, 0, 0), v=(3000, 3000, 0), s=(0, 1), time=1000., prob=10.)
            continue

        s.process(n)

        r = s.events
        return

    def test4(self):
        'Neutrons passed the surface'

        from sample_component import Sample

        rq = lambda Q: Q * 0.5 + 1

        s = Sample('sample', 5, 5, rq)

        from mcni import neutron_buffer, neutron
        N = 8
        n = neutron_buffer(N)
        for i in range(N):
            n[i] = neutron(r=(0, 0, 8), v=(0, 3000, 30000), s=(0, 1), time=1000., prob=10.)
            continue

        s.process(n)

        r = s.events
        return
        

    pass # end of TestCase


def pysuite():
    suite1 = unittest.makeSuite(TestCase)
    return unittest.TestSuite( (suite1,) )


def main():
    #debug.activate()
    pytests = pysuite()
    alltests = unittest.TestSuite( (pytests, ) )
    res = unittest.TextTestRunner(verbosity=2).run(alltests)
    import sys; sys.exit(not res.wasSuccessful())

    return


if __name__ == "__main__": 
    interactive = True
    main()

    
# version
__id__ = "$Id$"

# End of file 
