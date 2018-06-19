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
            print self.c
        else:
            self.c  = initial_covariance

        if accel_bias is not None:
            self.c[2,2] = accel_bias

    def bayard(self, t):
        k1 = 0.05 # 6/5!
        k2 = 0.25 # 6/4!
        
        c00 = k1*self.q2 * t**5 + k2*self.c[2,2] * t**4 + (self.c[1,2] + self.q1/6.0)*t**3 + (self.c[0,2] + self.c[1,1]) * t**2 + (self.c[0,1] + self.q0) * t + self.c[0,0]

        k3 = 0.5
        #c11 = k2*self.q2 * t**4 + self.c[2,2] * t**3 + (3.0 * self.c[1,2] + self.q1) * t**2 + (self.c[0,2] + self.c[1,1]) * t + 2.0 * self.c[0,1] + self.q0
        c11 = self.q2 * t**3 + self.c[2,2] * t**2 + (2*self.c[1,2] + self.q1) * t + self.c[1,1] #  + self.b**2

        return c00, c11


if __name__ == '__main__':
    honeywell_mimu_qa2000 = Accelerometer(sampling_frequency    = 200.0, # Hz
                                          velocity_read_noise   = 0.003/3.0, # m/s
                                          velocity_random_walk  = 0.0,
                                          accel_random_walk     = (181.0/3.0 * 9.81e-6)**2 / (365*24*3600.0), # ug over 12 months
                                          accel_bias            = 0.0027 / 60.0) # m/s/count to m/s2 (velocity quantization)
    lsm6dsl1              = Accelerometer(sampling_frequency    = 1666.0,
                                          position_random_walk  = 0.0,
                                          accel_random_walk     = (80 * 9.81e-6)**2, # ug/sqrt(Hz) to m2/s5
                                          accel_meas_variance   = 1800.0 * 9.81e-6, # ug(RMS) to m/s2
                                          accel_bias            = 40.0 * 9.81e-3)
    lsm6dsl2              = Accelerometer(sampling_frequency    = 1666.0,
                                          position_random_walk  = 0.0,
                                          accel_random_walk     = (80 * 9.81e-6)**2, # ug/sqrt(Hz) to m2/s5
                                          accel_meas_variance   = 2000.0 * 9.81e-6, # ug(RMS) to m/s2
                                          accel_bias            = 40.0 * 9.81e-3)
    lsm6dsl3              = Accelerometer(sampling_frequency    = 1666.0,
                                          position_random_walk  = 0.0,
                                          accel_random_walk     = (90 * 9.81e-6)**2, # ug/sqrt(Hz) to m2/s5
                                          accel_meas_variance   = 2400.0 * 9.81e-6, # ug(RMS) to m/s2
                                          accel_bias            = 40.0 * 9.81e-3)
    lsm6dsl4              = Accelerometer(sampling_frequency    = 1666.0,
                                          position_random_walk  = 0.0,
                                          accel_random_walk     = (130 * 9.81e-6)**2, # ug/sqrt(Hz) to m2/s5
                                          accel_meas_variance   = 3000.0 * 9.81e-6, # ug(RMS) to m/s2
                                          accel_bias            = 40.0 * 9.81e-3)
    vib                   = Accelerometer(sampling_frequency    = 1300.0,
                                          accel_random_walk     = (2.7 * 9.8e-3)**2, # mg/sqrt(Hz) to m2/s5
                                          accel_meas_variance   = 41.0, # m/s2
                                          accel_bias            = 1.5 * 9.81e-3,
                                          position_random_walk  = 0.0)

    time = np.arange(0, 60, 0.001)
    mimu_pos_cov  = np.zeros_like(time)
    mimu_vel_cov  = np.zeros_like(time)
    mode1_pos_cov = np.zeros_like(time)
    mode1_vel_cov = np.zeros_like(time)
    mode2_pos_cov = np.zeros_like(time)
    mode2_vel_cov = np.zeros_like(time)
    mode3_pos_cov = np.zeros_like(time)
    mode3_vel_cov = np.zeros_like(time)
    mode4_pos_cov = np.zeros_like(time)
    mode4_vel_cov = np.zeros_like(time)
    vib_pos_cov = np.zeros_like(time)
    vib_vel_cov = np.zeros_like(time)    
    mimu_pos_cov[0] = honeywell_mimu_qa2000.c[0,0]
    mimu_vel_cov[0] = honeywell_mimu_qa2000.c[1,1]
    mode1_pos_cov[0] = lsm6dsl1.c[0,0]
    mode1_vel_cov[0] = lsm6dsl1.c[1,1]
    mode2_pos_cov[0] = lsm6dsl2.c[0,0]
    mode2_vel_cov[0] = lsm6dsl2.c[1,1]
    mode3_pos_cov[0] = lsm6dsl3.c[0,0]
    mode3_vel_cov[0] = lsm6dsl3.c[1,1]
    mode4_pos_cov[0] = lsm6dsl4.c[0,0]
    mode4_vel_cov[0] = lsm6dsl4.c[1,1]
    vib_pos_cov[0] = vib.c[0,0]
    vib_vel_cov[0] = vib.c[1,1]         

    ii = 1
    for t in time[1:]:
        mimu_pos_cov[ii], mimu_vel_cov[ii]   = honeywell_mimu_qa2000.bayard(t)
        mode1_pos_cov[ii], mode1_vel_cov[ii] = lsm6dsl1.bayard(t)
        mode2_pos_cov[ii], mode2_vel_cov[ii] = lsm6dsl2.bayard(t)
        mode3_pos_cov[ii], mode3_vel_cov[ii] = lsm6dsl3.bayard(t)
        mode4_pos_cov[ii], mode4_vel_cov[ii] = lsm6dsl4.bayard(t)
        vib_pos_cov[ii], vib_vel_cov[ii] = vib.bayard(t)
        ii += 1
        
    fig = plt.figure()
    fig.suptitle("INS position and velocity error growth")
    ax1 = fig.add_subplot(211)
    #ax1.set_title("INS position error growth")
    ax1.plot(time, np.sqrt(mimu_pos_cov), label="Honeywell QA-2000")
    ax1.plot(time, np.sqrt(mode1_pos_cov), label="LSM6DSL Mode 1")
    ax1.plot(time, np.sqrt(mode2_pos_cov), label="LSM6DSL Mode 2")
    ax1.plot(time, np.sqrt(mode3_pos_cov), label="LSM6DSL Mode 3")
    ax1.plot(time, np.sqrt(mode4_pos_cov), label="LSM6DSL Mode 4")
    ax1.plot(time, np.sqrt(vib_pos_cov),   label="ADXL377 200g vibrometer")
    #ax1.set_xlabel("Time (seconds)")
    ax1.set_ylabel("m")
    ax1.set_yscale("log")
    ax1.set_xscale("log", nonposx='clip')
    ax1.legend()
    ax1.grid(True)

    ax2 = fig.add_subplot(212, sharex=ax1)
    #ax2.set_title("INS velocity error growth")
    ax2.plot(time, np.sqrt(mimu_vel_cov), label="Honeywell QA-2000")
    ax2.plot(time, np.sqrt(mode1_vel_cov), label="LSM6DSL Mode 1")
    ax2.plot(time, np.sqrt(mode2_vel_cov), label="LSM6DSL Mode 2")
    ax2.plot(time, np.sqrt(mode3_vel_cov), label="LSM6DSL Mode 3")
    ax2.plot(time, np.sqrt(mode4_vel_cov), label="LSM6DSL Mode 4")
    ax2.plot(time, np.sqrt(vib_vel_cov),   label="ADXL377 200g vibrometer")    
    ax2.set_xlabel("Time (seconds)")
    ax2.set_ylabel("m/s")
    ax2.set_yscale("log")
    ax2.set_xscale("log", nonposx='clip')
    ax2.set_xlim((1e-3, 1e2))
    ax2.grid(True)
    
    plt.show()
