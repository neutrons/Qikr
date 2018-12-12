import numpy as np
import mcvine, mcvine.components
from mcni.utils import conversion

def qikr(wl_0=1., center_wl=13.05):
    """
        @param wl_min: wavelength at the lower end of the bandwidth
    """
    running_length = 0
    k0 = 2*np.pi/wl_0
    e0_max = conversion.k2e(k0)
    kmin = 2*np.pi/40.0
    e0_min = conversion.k2e(kmin)
    instrument = mcvine.instrument()

    # Source - includes straight tube
    beam_width = 0.0266
    tube_length = 1.0
    mod = mcvine.components.sources.SNS_source_r1('moderator',
                                                  S_filename='./SNS-STS-ROT1-15_sp.dat',
                                                  width=beam_width, height=0.0266, dist=tube_length,
                                                  xw=beam_width, yh=0.0327, Emin=0, Emax=e0_max)
    instrument.append(mod, position=(0,0,0))
    #running_length += 0.001

    # Source Tube  ##################################################################
    #tube = mcvine.components.optics.Guide_channeled(name='source_tube',
    #                                                w1=beam_width, h1=0.0266, w2=beam_width, h2=0.0327,
    #                                                l=tube_length, mx=0, my=0)
    running_length += tube_length + 0.0001

    # Ballistic Guide 1 #############################################################
    # From 100 cm to 620 cm
    guide_1_length = 5.2
    guide1 = mcvine.components.optics.Guide_channeled(name='guide_1',
                                                      w1=beam_width, h1=0.0327, w2=beam_width, h2=0.0642,
                                                      l=guide_1_length, mx=5, my=5)
    instrument.append(guide1, position=(0,0,running_length))
    print ("Ballistic guide 1: %s m" % running_length)
    running_length += guide_1_length + 0.0001

    # T0 chopper ###################################################################
    # From 620 cm to 664.8 cm
    #t0_chopper = mcvine.components.optics.Vertical_T0(name='t0chopper', len=0.448,
    #                                                  w1=beam_width, w2=beam_width, nu=120, delta=0.0, tc=0., ymin=-.045, ymax=0.045)
    #instrument.append(t0_chopper, position=(0,0,6.2))
    running_length += 0.448

    # Ballistic Guide 2 #############################################################
    # From 664.8 cm to 714.8 cm
    guide_2_length = 0.5
    guide2 = mcvine.components.optics.Guide_channeled(name='guide_2',
                                                      w1=beam_width, h1=0.067, w2=beam_width, h2=0.07,
                                                      l=guide_2_length, mx=5, my=5)
    instrument.append(guide2, position=(0,0,running_length))
    print ("Ballistic guide 2: %s m" % running_length)
    running_length += guide_2_length + 0.0001

    # Bandwidth chopper ############################################################
    # From 714.8 cm to 722.1 cm
    instrument.append(mcvine.components.monitors.NeutronToStorage('save_pre', 'beam_pre_chopper.neutrons'),
                      position=(0,0,running_length))
    running_length += 0.0001
    b_chopper_length = 0.073
    b_chopper_dist = running_length + b_chopper_length/2.0  # center of chopper
    print ("Bandwidth chopper: %s m" % b_chopper_dist)

    aperture_bottom = 0.2689
    aperture_top = 0.33949
    radius_to_beam = 0.30390
    open_angle = 111.964

    kmid = 2*np.pi/center_wl
    v = conversion.k2v(kmid)
    delay = b_chopper_dist/v
    disk_chopper = mcvine.components.optics.DiskChopper_v2(name='band_chopper', 
                                                           yheight=aperture_top-aperture_bottom,
                                                           nslit=1,
                                                           radius=aperture_top,
                                                           theta_0=open_angle,
                                                           delay=delay, nu=7.5)
    running_length += 0.0001
    instrument.append(disk_chopper, position=(0,0, b_chopper_dist))
    instrument.append(mcvine.components.monitors.NeutronToStorage('save_post', 'beam_post_chopper.neutrons'),
                      position=(0,0,running_length))
    running_length += b_chopper_length + 0.0001

    # Ballistic Guide 3 #############################################################
    # From 722.1 cm to 1264 cm
    guide_3_length = 5.419
    guide3 = mcvine.components.optics.Guide_channeled(name='guide_3',
                                                      w1=beam_width, h1=0.07, w2=beam_width, h2=0.1033,
                                                      l=guide_3_length, mx=5, my=5)
    instrument.append(guide3, position=(0,0,running_length))
    running_length += guide_3_length + 0.0001


    # Tapered guide #################################################################
    # From 1264 cm to 1599 cm
    guide_4_length = 3.35
    guide4 = mcvine.components.optics.Guide_channeled(name='guide_4',
                                                      w1=beam_width, h1=0.1033, w2=beam_width, h2=0.02,
                                                      l=guide_4_length, mx=5, my=5)
    instrument.append(guide4, position=(0,0,running_length))
    running_length += guide_4_length + 0.0001

    instrument.append(mcvine.components.monitors.NeutronToStorage('save', 'beam.neutrons'),
                      position=(0,0,running_length))
    running_length += 0.0001
    print("End of instrument: %s m" % running_length)

    return instrument
