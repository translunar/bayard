"""Check Gyroscope.bayard() against a direct transcription of the reference
MATLAB in references/bayard_calc.m.

The gyro and star tracker values are the `jpl_mimu` and `st_bct` cases from
references/bayard_method.m.
"""

import numpy as np

from bayard import Gyroscope


D2R  = np.pi / 180.0
AS2D = 1.0 / 3600.0

# references/bayard_method.m, gyro.jpl_mimu
JPL_MIMU_BIAS_STABILITY = 0.05 / 3 * AS2D * D2R   # rad/s
JPL_MIMU_RANDOM_WALK    = 0.025 / 3 * 1 / 60 * D2R  # rad/sqrt(s)

# references/bayard_method.m, st_bct (BCT Nano Star Tracker)
ST_NEA   = 333e-6            # rad, 1-sigma
ST_DELTA = 0.2               # s
ST_R     = ST_DELTA * ST_NEA**2   # rad2 s
ST_B     = 60 * AS2D * D2R   # rad, 1-sigma

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


def make_gyro():
    return Gyroscope(angle_random_walk             = JPL_MIMU_RANDOM_WALK**2,
                     bias_stability                = JPL_MIMU_BIAS_STABILITY**2 / 3600.0,
                     attitude_meas_sigma           = ST_NEA,
                     attitude_meas_sampling_period = ST_DELTA,
                     attitude_meas_bias            = ST_B)


def test_initial_covariance_matches_matlab():
    """c[0,0], c[0,1] and c[1,1] must equal p11, p12 and p22."""
    g = make_gyro()

    q1 = JPL_MIMU_RANDOM_WALK**2
    q2 = JPL_MIMU_BIAS_STABILITY**2 / 3600.0
    l  = np.sqrt(q1 + 2 * np.sqrt(ST_R * q2))

    assert np.isclose(g.c[0, 0], np.sqrt(ST_R) * l, rtol=1e-12)   # p11
    assert np.isclose(g.c[0, 1], np.sqrt(ST_R * q2), rtol=1e-12)  # p12
    assert np.isclose(g.c[1, 1], np.sqrt(q2) * l, rtol=1e-12)     # p22
    assert g.c[1, 0] == g.c[0, 1]


def test_bayard_matches_matlab():
    """p(t) must match bayard_calc.m at t = 1, 60 and 3600 s."""
    g = make_gyro()

    for t in TIMES:
        expected = bayard_calc(JPL_MIMU_RANDOM_WALK, JPL_MIMU_BIAS_STABILITY,
                               ST_R, ST_B, t)
        assert np.isclose(g.bayard(t), expected, rtol=1e-12), \
            "t = %g: got %r, expected %r" % (t, g.bayard(t), expected)


def test_bayard_is_monotonic_and_starts_at_p11():
    g = make_gyro()

    # p(0) is p11 + b^2
    assert np.isclose(g.bayard(0.0), g.c[0, 0] + ST_B**2, rtol=1e-12)

    values = [g.bayard(t) for t in TIMES]
    assert all(a < b for a, b in zip(values, values[1:]))


if __name__ == '__main__':
    test_initial_covariance_matches_matlab()
    test_bayard_matches_matlab()
    test_bayard_is_monotonic_and_starts_at_p11()
    print("all tests passed")
