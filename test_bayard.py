"""Check bayard.py against the reference implementations.

Gyroscope.bayard() is checked against a direct transcription of the MATLAB in
references/bayard_calc.m, using the `jpl_mimu` and `st_bct` values from
references/bayard_method.m.

Accelerometer.bayard() has no MATLAB to check against -- bayard_calc.m and
bayard_method.m cover the gyro only -- so it is checked against a direct
transcription of equation (1.6) of references/bayard2000.pdf.
"""

import unittest
from math import factorial

import numpy as np

from bayard import Accelerometer, Gyroscope


D2R  = np.pi / 180.0
AS2D = 1.0 / 3600.0

# references/bayard_method.m, gyro.jpl_mimu
JPL_MIMU_BIAS_STABILITY = 0.05 / 3 * AS2D * D2R     # rad/s
JPL_MIMU_RANDOM_WALK    = 0.025 / 3 * 1 / 60 * D2R  # rad/sqrt(s)

# references/bayard_method.m, st_bct (BCT Nano Star Tracker)
ST_NEA   = 333e-6                 # rad, 1-sigma
ST_DELTA = 0.2                    # s
ST_R     = ST_DELTA * ST_NEA**2   # rad2 s
ST_B     = 60 * AS2D * D2R        # rad, 1-sigma

TIMES = (1.0, 60.0, 3600.0)


def bayard_calc(random_walk, bias_stability, r, b, dt):
    """Direct transcription of references/bayard_calc.m (returns p, not sqrt(p))."""
    q1 = random_walk**2
    q2 = bias_stability**2 / 3600.0

    l = np.sqrt(q1 + 2 * np.sqrt(r * q2))

    p11 = np.sqrt(r) * l
    p12 = np.sqrt(r * q2)
    p22 = np.sqrt(q2) * l

    return q2 / 3 * dt**3 + p22 * dt**2 + (2 * p12 + q1) * dt + p11 + b**2


class TestGyroscope(unittest.TestCase):
    """Gyroscope must reproduce bayard_calc.m for the jpl_mimu / st_bct case."""

    RTOL = 1e-12

    def setUp(self):
        self.gyro = Gyroscope(angle_random_walk             = JPL_MIMU_RANDOM_WALK**2,
                              bias_stability                = JPL_MIMU_BIAS_STABILITY**2 / 3600.0,
                              attitude_meas_sigma           = ST_NEA,
                              attitude_meas_sampling_period = ST_DELTA,
                              attitude_meas_bias            = ST_B)

        self.q1 = JPL_MIMU_RANDOM_WALK**2                    # q1 in bayard_calc.m
        self.q2 = JPL_MIMU_BIAS_STABILITY**2 / 3600.0        # q2 in bayard_calc.m
        self.l  = np.sqrt(self.q1 + 2 * np.sqrt(ST_R * self.q2))

    def assertClose(self, actual, expected, msg = None):
        """Purely relative comparison.

        atol must be 0.0. These covariances run down to 1e-15 rad2, so
        numpy's default atol of 1e-08 would swamp them entirely and the
        comparison would pass against anything, zero included.
        """
        self.assertTrue(np.isclose(actual, expected, rtol = self.RTOL, atol = 0.0),
                        msg or "got %r, expected %r" % (actual, expected))

    def test_assertclose_rejects_zero(self):
        """Guard against reintroducing an absolute tolerance.

        c[1,1] is around 1e-15, so with numpy's default atol of 1e-08 this
        comparison against zero would succeed and every other assertion in
        this file would be vacuous.
        """
        with self.assertRaises(AssertionError):
            self.assertClose(self.gyro.c[1, 1], 0.0)

    def test_initial_covariance_matches_matlab(self):
        """c[0,0], c[0,1] and c[1,1] must equal p11, p12 and p22."""
        self.assertClose(self.gyro.c[0, 0], np.sqrt(ST_R) * self.l)     # p11
        self.assertClose(self.gyro.c[0, 1], np.sqrt(ST_R * self.q2))    # p12
        self.assertClose(self.gyro.c[1, 1], np.sqrt(self.q2) * self.l)  # p22

    def test_initial_covariance_is_symmetric_two_by_two(self):
        self.assertEqual(self.gyro.c.shape, (2, 2))
        self.assertEqual(self.gyro.c[1, 0], self.gyro.c[0, 1])

    def test_bayard_matches_matlab(self):
        """p(t) must match bayard_calc.m at t = 1, 60 and 3600 s."""
        for t in TIMES:
            expected = bayard_calc(JPL_MIMU_RANDOM_WALK, JPL_MIMU_BIAS_STABILITY,
                                   ST_R, ST_B, t)
            self.assertClose(self.gyro.bayard(t), expected,
                             "t = %g: got %r, expected %r" % (t, self.gyro.bayard(t), expected))

    def test_bayard_starts_at_p11_plus_b_squared(self):
        self.assertClose(self.gyro.bayard(0.0), self.gyro.c[0, 0] + ST_B**2)

    def test_bayard_is_monotonic(self):
        values = [self.gyro.bayard(t) for t in TIMES]
        for earlier, later in zip(values, values[1:]):
            self.assertLess(earlier, later)


def bayard2000_eq_1_6(q0, q1, q2, c, t):
    """Direct transcription of equation (1.6) of references/bayard2000.pdf.

    The memo indexes the position/rate/accel error vector from 1, so its
    c11/c12/c13/c22/c23/c33 are c[0,0]/c[0,1]/c[0,2]/c[1,1]/c[1,2]/c[2,2] here.
    Written with the factorials left in, to stay close to the printed equation.
    """
    c11_0, c12_0, c13_0 = c[0, 0], c[0, 1], c[0, 2]
    c22_0, c23_0, c33_0 = c[1, 1], c[1, 2], c[2, 2]

    return (6 * q2 / factorial(5)) * t**5 \
        + 6 * (c33_0 / factorial(4)) * t**4 \
        + ((6 * c23_0 + 2 * q1) / factorial(3)) * t**3 \
        + 2 * ((c13_0 + c22_0) / factorial(2)) * t**2 \
        + (2 * c12_0 + q0) * t \
        + c11_0


def velocity_propagation(q1, q2, c, t):
    """P11(t) for the triple integrator.

    The memo's section 1.2 promises a velocity covariance c22(t) but never
    gives one, so this is derived rather than transcribed:
    P(t) = Phi P(0) Phi' + integral of Phi Q Phi', with Phi the triple
    integrator state transition matrix [[1,t,t2/2],[0,1,t],[0,0,1]].
    """
    return c[1, 1] + (2 * c[1, 2] + q1) * t + c[2, 2] * t**2 + q2 * t**3 / 3.0


class TestAccelerometer(unittest.TestCase):
    """Accelerometer must reproduce equation (1.6) of the Bayard memo."""

    RTOL = 1e-12

    # Arbitrary but fully populated: every initial covariance term must be
    # nonzero, or a dropped term goes unnoticed. In particular c[0,1] must be
    # nonzero to pin down the factor of 2 on the linear term.
    Q0 = 3.0e-7   # position random walk, m2/s
    Q1 = 5.0e-9   # velocity random walk, m2/s3
    Q2 = 7.0e-13  # accel random walk,    m2/s5

    C0 = np.array([[4.0e-4, 3.0e-5, 2.0e-7],
                   [3.0e-5, 1.2e-4, 5.0e-7],
                   [2.0e-7, 5.0e-7, 9.0e-9]])

    def setUp(self):
        self.accel = Accelerometer(position_random_walk = self.Q0,
                                   velocity_random_walk = self.Q1,
                                   accel_random_walk    = self.Q2,
                                   initial_covariance   = self.C0.copy())

    def assertClose(self, actual, expected, msg = None):
        """Purely relative; see TestGyroscope.assertClose for why atol is 0.0."""
        self.assertTrue(np.isclose(actual, expected, rtol = self.RTOL, atol = 0.0),
                        msg or "got %r, expected %r" % (actual, expected))

    def test_position_matches_memo_eq_1_6(self):
        for t in TIMES:
            expected = bayard2000_eq_1_6(self.Q0, self.Q1, self.Q2, self.C0, t)
            c00, _ = self.accel.bayard(t)
            self.assertClose(c00, expected,
                             "t = %g: got %r, expected %r" % (t, c00, expected))

    def test_position_linear_term_carries_factor_of_two(self):
        """The t coefficient is (2 c12(0) + q0), not (c12(0) + q0)."""
        # Isolate the linear term: with only c[0,1] and q0 nonzero, and the
        # higher-order initial terms zeroed, c00(t) - c00(0) is exactly
        # (2 c12 + q0) t.
        c = np.zeros((3, 3))
        c[0, 1] = c[1, 0] = 3.0e-5
        accel = Accelerometer(position_random_walk = self.Q0,
                              velocity_random_walk = 0.0,
                              accel_random_walk    = 0.0,
                              initial_covariance   = c)
        t = 10.0
        c00, _ = accel.bayard(t)
        self.assertClose(c00 - c[0, 0], (2 * c[0, 1] + self.Q0) * t)

    def test_velocity_matches_derivation(self):
        for t in TIMES:
            expected = velocity_propagation(self.Q1, self.Q2, self.C0, t)
            _, c11 = self.accel.bayard(t)
            self.assertClose(c11, expected,
                             "t = %g: got %r, expected %r" % (t, c11, expected))

    def test_bayard_starts_at_initial_covariance(self):
        c00, c11 = self.accel.bayard(0.0)
        self.assertClose(c00, self.C0[0, 0])
        self.assertClose(c11, self.C0[1, 1])


if __name__ == '__main__':
    unittest.main()
