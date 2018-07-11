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


def convert_errors(values, variance = True, convert = True):
    """Convert variance values to standard deviations and potentially
    also to ft, ft/sec.
    """
    FEET_PER_METER = 3.28084
    
    if variance:
        if convert: return (np.sqrt(values) * FEET_PER_METER)**2
        else:       return values
    else:
        if convert: return np.sqrt(values) * FEET_PER_METER
        else:       return np.sqrt(values)


def ylabels(variance = True, convert = True):
    if variance:
        if convert: return ("ft2", "ft2/s2")
        else:       return ("m2",  "m2/s2")
    else:
        if convert: return ("ft",  "ft/s")
        else:       return ("m",   "m/s")
    


def plot_errors(time, pos_cov, vel_cov,
                which=('lsm6dsl1', 'lsm6dsl3', 'lsm6dsl4', 'vibxy', 'vibz', 'vibtemp'),
                labels = {'honeywell': "Honeywell",
                          'lsm6dsl1' : "Mode 1-2",
                          'lsm6dsl3' : "Mode 3",
                          'lsm6dsl4' : "Mode 4",
                          'vibxy'    : "Vibrometer xy",
                          'vibz'     : "Vibrometer z",
                          'vibtemp'  : "Vibrometer z + 30K"},
                colors = {'honeywell': 'red',
                          'lsm6dsl1' : 'blue',
                          'lsm6dsl3' : 'orange',
                          'lsm6dsl4' : 'green',
                          'vibxy'    : 'purple',
                          'vibz'     : 'indigo',
                          'vibtemp'  : 'brown'},
                filename = 'bayard_accel.pdf',
                variance = True,
                xscale   = None,
                yscale   = None,
                xlim     = None,
                pos_ylim = None,
                vel_ylim = None,
                convert  = True): # convert to ft, ft/sec

    
    
    pdf = PdfPages(filename)
    fig = plt.figure(figsize=(6,8))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)
    ax1.grid(True)
    ax2.grid(True)
    plt.setp(ax1.get_xticklabels(), visible=False)
    fig.subplots_adjust(hspace=0)

    ax2.set_xlabel("time (s)")

    yl = ylabels(variance, convert)
    ax1.set_ylabel(yl[0])
    ax2.set_ylabel(yl[1])

    for key in which:
        y1 = convert_errors(pos_cov[key], variance, convert)
        y2 = convert_errors(vel_cov[key], variance, convert)
        
        ax1.plot(time, y1, label=labels[key], color=colors[key])
        ax2.plot(time, y2, label=labels[key], color=colors[key])

    if xscale is not None: ax1.set_xscale(xscale)
    if yscale is not None:
        ax1.set_yscale(yscale)
        ax2.set_yscale(yscale)

    if pos_ylim: ax1.set_ylim(pos_ylim)
    if vel_ylim: ax2.set_ylim(vel_ylim)
    if xlim:
        ax1.set_xlim(xlim)
    else:
        ax1.set_xlim(xmax=time[-1])
        
    ax1.legend()
    pdf.savefig(fig)
    pdf.close()

    return fig
                            

class Gyroscope(object):
    def __init__(self,
                 angle_random_walk             = 5.8761e-12, # r2/s
                 bias_stability                = 1.8138e-18, # r2/s3
                 angle_variance                = None,
                 bias_variance                 = None,                 
                 attitude_meas_variance        = None,  # rad
                 attitude_meas_sampling_period = 1.0,   # s
                 attitude_meas_bias            = 0.0,   # rad
                 initial_covariance            = None): # rad
        self.q0 = angle_random_walk
        self.q1 = bias_stability

        if initial_covariance is not None:
            self.c = initial_covariance
        else:
            if attitude_meas_variance is not None:
                r = attitude_meas_sampling_period * attitude_meas_variance
                l = np.sqrt(self.q1 + 2.0 * np.sqrt(r * self.q1))
                self.c = np.zeros((3,3))
                self.c[0,0] = np.sqrt(r) * l
                self.c[0,1] = np.sqrt(r * self.q1)
                self.c[1,1] = np.sqrt(self.q1)
                self.c[1,0] = self.c[0,1]
            else:
                print("Warning: Cross-correlation for angle and bias is 0")
                self.c[0,1] = self[c,10] = 0.0

            if angle_variance is not None:
                self.c[0,0] = angle_variance
            if bias_variance is not None:
                self.c[1,1] = bias_variance

        if attitude_meas_bias is None:
            self.b = 0.0
        else:
            self.b = attitude_meas_bias

    def bayard(self, t):
        return self.q1 * t**3 / 3.0 + self.c[1,1]*t**2 + (2.0 * self.c[0,1] + self.q0)*t + self.c[0,0] + self.b
                 

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
                 accel_bias_variance         = None,
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

            print ("Warning: Cross-correlation for position and acceleration is 0")
            self.c[0,2] = self.c[2,0] = 0.0 # FIXME: not really sure what this should be
        else:
            self.c  = initial_covariance

        if accel_bias_variance is not None:
            self.c[2,2] = accel_bias_variance
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
                                          accel_bias_variance   = (0.0027 / 60.0)**2) # m/s/count to m2/s4 (velocity quantization)
    position_variance = honeywell_mimu_qa2000.c[0,0] # 0.00036 (matches 3D analysis)

    acc = {}
    acc['honeywell']      = honeywell_mimu_qa2000
    acc['lsm6dsl1']       = Accelerometer(velocity_variance     = velocity_variance,
                                          position_variance     = position_variance,
                                          position_random_walk  = 0.0,
                                          accel_random_walk     = (80 * 9.81e-6)**2, # ug/sqrt(Hz) to m2/s5
                                          #accel_meas_variance   = 1800.0 * 9.81e-6, # ug(RMS) to m/s2
                                          accel_bias_variance   = (40.0 * 9.81e-3)**2)
    acc['lsm6dsl2']       = Accelerometer(velocity_variance     = velocity_variance,
                                          position_variance     = position_variance,
                                          position_random_walk  = 0.0,
                                          accel_random_walk     = (80 * 9.81e-6)**2, # ug/sqrt(Hz) to m2/s5
                                          #accel_meas_variance   = 2000.0 * 9.81e-6, # ug(RMS) to m/s2
                                          accel_bias_variance   = (40.0 * 9.81e-3)**2)
    acc['lsm6dsl3']       = Accelerometer(velocity_variance     = velocity_variance,
                                          position_variance     = position_variance,
                                          position_random_walk  = 0.0,
                                          accel_random_walk     = (90 * 9.81e-6)**2, # ug/sqrt(Hz) to m2/s5
                                          #accel_meas_variance   = 2400.0 * 9.81e-6, # ug(RMS) to m/s2
                                          accel_bias_variance   = (40.0 * 9.81e-3)**2)
    acc['lsm6dsl4']       = Accelerometer(velocity_variance     = velocity_variance,
                                          position_variance     = position_variance,
                                          position_random_walk  = 0.0,
                                          accel_random_walk     = (130 * 9.81e-6)**2, # ug/sqrt(Hz) to m2/s5
                                          #accel_meas_variance   = 3000.0 * 9.81e-6, # ug(RMS) to m/s2
                                          accel_bias_variance  = (40.0 * 9.81e-3)**2)
    acc['vibxy']          = Accelerometer(velocity_variance     = velocity_variance,
                                          position_variance     = position_variance,
                                          accel_random_walk     = (2.7 * 9.8e-3)**2, # mg/sqrt(Hz) to m2/s5
                                          #accel_meas_variance   = 41.0, # m/s2
                                          accel_bias_variance   = (12 * 9.81e-3)**2,
                                          position_random_walk  = 0.0)
    acc['vibz']           = Accelerometer(velocity_variance     = velocity_variance,
                                          position_variance     = position_variance,
                                          accel_random_walk     = (4.3 * 9.8e-3)**2, # mg/sqrt(Hz) to m2/s5
                                          #accel_meas_variance   = 41.0, # m/s2
                                          accel_bias_variance   = (30 * 9.81e-3)**2,
                                          position_random_walk  = 0.0)
    acc['vibtemp']        = Accelerometer(velocity_variance     = velocity_variance,
                                          position_variance     = position_variance,
                                          accel_random_walk     = (4.3 * 9.8e-3)**2, # mg/sqrt(Hz) to m2/s5
                                          #accel_meas_variance   = 41.0, # m/s2
                                          accel_bias_variance   = (900 * 9.81e-3)**2,
                                          position_random_walk  = 0.0)


    #import pdb
    #pdb.set_trace()
    gyro = {}
    gyro['lsm6dsl']   = Gyroscope(angle_random_walk = 0.004,
                                  bias_stability    = (3.0 * np.pi/180.0),
                                  attitude_meas_variance = (333e-6)**2,
                                  attitude_meas_bias     = (333e-6)**2,
                                  attitude_meas_sampling_period = 0.001)
    gyro['honeywell'] = Gyroscope(angle_random_walk = 5.8761e-12, # r2/s
                                  bias_stability    = 1.8138e-18, # r2/s3,
                                  attitude_meas_variance = (333e-6)**2, #r2
                                  attitude_meas_bias     = (333e-6)**2, #r2
                                  attitude_meas_sampling_period = 0.5) # s
                                  
                                
    
    time = np.arange(0, 2.0, 0.001)#100.0, 0.001)

    att_cov = {}
    pos_cov = {}
    vel_cov = {}
    for aa in acc:
        ac = acc[aa]
        pos_cov[aa] = np.zeros_like(time)
        vel_cov[aa] = np.zeros_like(time)
    for gg in gyro:
        gc = gyro[gg]
        att_cov[gg] = np.zeros_like(time)

        att_cov[gg][0] = gc.c[0,0]

        ii = 1
        for t in time[1:]:
            att_cov[gg][ii] = gc.bayard(t)
            ii += 1
    
    
    for aa in acc:
        ac = acc[aa]

        pos_cov[aa][0] = ac.c[0,0]
        vel_cov[aa][0] = ac.c[1,1]

        ii = 1
        for t in time[1:]:
            pos_cov[aa][ii], vel_cov[aa][ii] = ac.bayard(t)
            ii += 1

    fig = plot_errors(time, pos_cov, vel_cov, filename='bayard_accel_variance.pdf', variance = True, xscale = 'log', yscale = 'log', convert=True, xlim=(0.001, 100.0))
    fig = plot_errors(time, pos_cov, vel_cov, filename='bayard_accel_sigma_1s_no_vibtemp.pdf', variance = False, convert=True, which=('lsm6dsl1', 'lsm6dsl3', 'lsm6dsl4', 'vibxy', 'vibz'))
    fig = plot_errors(time, pos_cov, vel_cov, filename='bayard_accel_sigma_1s.pdf', variance = False, yscale = 'log', convert=True)
    #fig = plot_errors(time, pos_cov, vel_cov, filename='bayard_accel_sigma.pdf', variance = False, xscale = 'log', yscale = 'log', convert=True)
    #fig = plot_errors(time, pos_cov, vel_cov, filename='bayard_accel_sigma_1s.pdf', variance = False, xlim=(0.0, 1.0), xscale = 'log', yscale = 'log', vel_ylim=(0.03, 50.0), pos_ylim=(0.02, 25.0), convert=True)
    #fig = plot_errors(time, pos_cov, vel_cov, filename='bayard_accel_sigma_5s.pdf', variance = False, xlim=(0.0, 2.0), which=('lsm6dsl4',), labels={'lsm6dsl4': "LSM6DSL"}, pos_ylim=(0.0, 3.0), vel_ylim=(0.0, 3.0), convert=True)
    
    plt.show()

