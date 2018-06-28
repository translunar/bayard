#!/usr/bin/python
"""For understanding covariance growth due to gyroscope and
accelerometer when only occasional attitude and position/velocity
measurements are available.

Reference:

* Bayard, D. 2000 Nov 14. First look analysis of IMU propagation
  errors for Mars Aerocapture Demonstration. Engineering Memorandum
  EM-3455-00-005, Jet Propulsion Laboratory, Pasadena, CA, US.

"""

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class Accelerometer(object):
    def __init__(self,
                 sampling_frequency          = 200.0, # Hz
                 velocity_sampling_freq      = None,
                 position_sampling_freq      = None,
                 velocity_meas_variance      = 1.0,
                 position_meas_variance      = 1.0,
                 accel_meas_variance         = 0.0027 / 60.0,
                 velocity_read_noise         = 0.003/3.0, # m/s
                 position_random_walk        = None, # overridden by velocity_read_noise
                 velocity_random_walk        = 0.0,
                 accel_random_walk           = (181.0/3.0 * 9.81e-6)**2 / (365*24*3600.0),
                 accel_bias                  = None,
                 velocity_bias               = 0.0, # not sure why this would ever be nonzero
                 velocity_variance           = None,
                 position_variance           = None,
                 initial_covariance          = None):

        if position_random_walk is None:
            self.q0 = velocity_read_noise**2 / sampling_frequency
        else:
            self.q0 = position_random_walk # m2/s
        self.q1 = velocity_random_walk # m2/s3
        self.q2 = accel_random_walk # m2/s5

        if initial_covariance is None:
            v_delta = 1.0 if velocity_sampling_freq is None else 1.0 / velocity_sampling_freq
            p_delta = 1.0 if position_sampling_freq is None else 1.0 / position_sampling_freq

            self.c  = np.zeros((3,3))
            
            # velocity covariance terms
            # r is in m2/s2 * s = m2/s
            # sqrt(r * q2) is in m2/s3
            # l is in m/s^1.5
            r       = v_delta * velocity_meas_variance # equivalent to star tracker per-axis covariance, but for velocity
            l       = np.sqrt(self.q1 + 2 * np.sqrt(r * self.q2))

            self.c[2,2] = np.sqrt(self.q2) * l               # m/s2.5 * m/s1.5 = m2/s4
            self.c[1,2] = self.c[2,1] = np.sqrt(r * self.q2) # sqrt(m2/s * m2/s5) = m2/s3 
            self.c[1,1] = np.sqrt(r) * l                     # m/sqrt(s) * m/s^1.5 = m2/s2
            
            # position covariance terms
            # s is in m2 * s = m2s
            # m is in m/sqrt(s)
            s       = p_delta * position_meas_variance
            m       = np.sqrt(self.q0 + 2 * np.sqrt(s * self.q1)) # m/sqrt(s)
            
            self.c[0,0] = np.sqrt(s) * m # m*sqrt(s) * m/sqrt(s) = m2
            self.c[0,1] = self.c[1,0] = np.sqrt(s * self.q1) # sqrt(m2s * m2/s3) = m2/s
            self.c[1,1] = np.sqrt(self.q1) * m # m/s1.5 * m/sqrt(s) = m2/s2

            self.c[0,2] = self.c[2,0] = 0.0 # FIXME: not really sure what this should be
        else:
            self.c  = initial_covariance

        if accel_bias is not None:
            self.c[2,2] = accel_bias**2
        if velocity_variance is not None:
            self.c[1,1] = velocity_variance
        if position_variance is not None:
            self.c[0,0] = position_variance

        print self.c

    def bayard(self, t):
        k1 = 0.05 # 6/5!
        k2 = 0.25 # 6/4!
        
        c00 = k1*self.q2 * t**5 + k2*self.c[2,2] * t**4 + (self.c[1,2] + self.q1/3.0)*t**3 + (self.c[0,2] + self.c[1,1]) * t**2 + (self.c[0,1] + self.q0) * t + self.c[0,0]

        k3 = 0.5
        #c11 = k2*self.q2 * t**4 + self.c[2,2] * t**3 + (3.0 * self.c[1,2] + self.q1) * t**2 + (self.c[0,2] + self.c[1,1]) * t + 2.0 * self.c[0,1] + self.q0
        c11 = self.q2/3.0 * t**3 + self.c[2,2] * t**2 + (2.0*self.c[1,2] + self.q1) * t + self.c[1,1] #  + self.b**2

        return c00, c11


if __name__ == '__main__':
    velocity_variance = 0.011**2 #3.6e-10 (matches 3D analysis)
    honeywell_mimu_qa2000 = Accelerometer(sampling_frequency    = 200.0, # Hz
                                          velocity_variance     = velocity_variance,
                                          velocity_read_noise   = 0.003/3.0, # m/s
                                          velocity_random_walk  = 0.0,
                                          accel_random_walk     = (181.0/3.0 * 9.81e-6)**2 / (365*24*3600.0), # ug over 12 months
                                          accel_bias            = 0.0027 / 60.0) # m/s/count to m/s2 (velocity quantization)
    position_variance = honeywell_mimu_qa2000.c[0,0] # 0.00036 (matches 3D analysis)

    acc = {}
    acc['honeywell']      = honeywell_mimu_qa2000
    acc['lsm6dsl1']       = Accelerometer(velocity_variance     = velocity_variance,
                                          position_variance     = position_variance,
                                          position_random_walk  = 0.0,
                                          accel_random_walk     = (80 * 9.81e-6)**2, # ug/sqrt(Hz) to m2/s5
                                          #accel_meas_variance   = 1800.0 * 9.81e-6, # ug(RMS) to m/s2
                                          accel_bias            = 40.0 * 9.81e-3)
    acc['lsm6dsl2']       = Accelerometer(velocity_variance     = velocity_variance,
                                          position_variance     = position_variance,
                                          position_random_walk  = 0.0,
                                          accel_random_walk     = (80 * 9.81e-6)**2, # ug/sqrt(Hz) to m2/s5
                                          #accel_meas_variance   = 2000.0 * 9.81e-6, # ug(RMS) to m/s2
                                          accel_bias            = 40.0 * 9.81e-3)
    acc['lsm6dsl3']       = Accelerometer(velocity_variance     = velocity_variance,
                                          position_variance     = position_variance,
                                          position_random_walk  = 0.0,
                                          accel_random_walk     = (90 * 9.81e-6)**2, # ug/sqrt(Hz) to m2/s5
                                          #accel_meas_variance   = 2400.0 * 9.81e-6, # ug(RMS) to m/s2
                                          accel_bias            = 40.0 * 9.81e-3)
    acc['lsm6dsl4']       = Accelerometer(velocity_variance     = velocity_variance,
                                          position_variance     = position_variance,
                                          position_random_walk  = 0.0,
                                          accel_random_walk     = (130 * 9.81e-6)**2, # ug/sqrt(Hz) to m2/s5
                                          #accel_meas_variance   = 3000.0 * 9.81e-6, # ug(RMS) to m/s2
                                          accel_bias            = 40.0 * 9.81e-3)
    acc['vibxy']          = Accelerometer(velocity_variance     = velocity_variance,
                                          position_variance     = position_variance,
                                          accel_random_walk     = (2.7 * 9.8e-3)**2, # mg/sqrt(Hz) to m2/s5
                                          #accel_meas_variance   = 41.0, # m/s2
                                          accel_bias            = 12 * 9.81e-3,
                                          position_random_walk  = 0.0)
    acc['vibz']           = Accelerometer(velocity_variance     = velocity_variance,
                                          position_variance     = position_variance,
                                          accel_random_walk     = (4.3 * 9.8e-3)**2, # mg/sqrt(Hz) to m2/s5
                                          #accel_meas_variance   = 41.0, # m/s2
                                          accel_bias            = 30 * 9.81e-3,
                                          position_random_walk  = 0.0)
    acc['vibtemp']        = Accelerometer(velocity_variance     = velocity_variance,
                                          position_variance     = position_variance,
                                          accel_random_walk     = (4.3 * 9.8e-3)**2, # mg/sqrt(Hz) to m2/s5
                                          #accel_meas_variance   = 41.0, # m/s2
                                          accel_bias            = 900 * 9.81e-3,
                                          position_random_walk  = 0.0)

    time = np.arange(0, 100.0, 0.001)

    pos_cov = {}
    vel_cov = {}
    for aa in acc:
        ac = acc[aa]
        pos_cov[aa] = np.zeros_like(time)
        vel_cov[aa] = np.zeros_like(time)
    
    
    
    for aa in acc:
        ac = acc[aa]

        pos_cov[aa][0] = ac.c[0,0]
        vel_cov[aa][0] = ac.c[1,1]

        ii = 1
        for t in time[1:]:
            pos_cov[aa][ii], vel_cov[aa][ii] = ac.bayard(t)
            ii += 1

    pdf = PdfPages('bayard.pdf')
    fig = plt.figure(figsize=(6,8))
    #fig.suptitle("INS position and velocity variance growth")
    ax1 = fig.add_subplot(211)
    #ax1.set_title("INS position error growth")
    ax1.plot(time, np.sqrt(pos_cov['honeywell']), label="Honeywell")
    ax1.plot(time, np.sqrt(pos_cov['lsm6dsl1']), label="Mode 1-2")
    ax1.plot(time, np.sqrt(pos_cov['lsm6dsl3']), label="Mode 3")
    ax1.plot(time, np.sqrt(pos_cov['lsm6dsl4']), label="Mode 4")
    ax1.plot(time, np.sqrt(pos_cov['vibxy']), label="Vibrometer xy")
    ax1.plot(time, np.sqrt(pos_cov['vibz']),  label="Vibrometer z")
    ax1.plot(time, np.sqrt(pos_cov['vibtemp']), label="Vibrometer z + 30K")
    #ax1.set_xlabel("Time (seconds)")
    ax1.set_ylabel("m2")
    ax1.set_yscale("log")
    ax1.set_xscale("log", nonposx='clip')
    ax1.grid(True)
    plt.setp(ax1.get_xticklabels(), visible=False)

    ax2 = fig.add_subplot(212, sharex=ax1)
    #ax2.set_title("INS velocity error growth")
    ax2.plot(time, np.sqrt(vel_cov['honeywell']), label="Honeywell")
    ax2.plot(time, np.sqrt(vel_cov['lsm6dsl1']), label="Mode 1-2")
    ax2.plot(time, np.sqrt(vel_cov['lsm6dsl3']), label="Mode 3")
    ax2.plot(time, np.sqrt(vel_cov['lsm6dsl4']), label="Mode 4")
    ax2.plot(time, np.sqrt(vel_cov['vibxy']), label="Vibrometer xy")
    ax2.plot(time, np.sqrt(vel_cov['vibz']),  label="Vibrometer z")
    ax2.plot(time, np.sqrt(vel_cov['vibtemp']), label="Vibrometer z + 30K")
    ax2.set_xlabel("Time (seconds)")
    ax2.set_ylabel("m2/s2")
    ax2.set_yscale("log")
    ax2.set_xscale("log", nonposx='clip')
    ax2.set_xlim((1e-3, 1e2))
    ax2.grid(True)
    ax2.legend()
    
    #plt.show()
    fig.subplots_adjust(hspace=0)
    pdf.savefig(fig)
    pdf.close()
    plt.show()

