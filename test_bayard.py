"""Check Gyroscope.bayard() against a direct transcription of the reference
MATLAB in references/bayard_calc.m.

The gyro and star tracker values are the `jpl_mimu` and `st_bct` cases from
references/bayard_method.m.
"""

import unittest

import numpy as np

from bayard import Gyroscope


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
        """Relative comparison; these covariances are far too small for assertAlmostEqual."""
        self.assertTrue(np.isclose(actual, expected, rtol = self.RTOL),
                        msg or "got %r, expected %r" % (actual, expected))

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


if __name__ == '__main__':
    unittest.main()
