"""Ice saturation, inverse and humidity checks independent of the liquid branch."""
import numpy as np
import pytest

from wsp2p import (
    T_from_e_ice, coeffs_ice, dln_esat_ice_dT, esat_ice_hpa,
    esat_water_hpa, frostpoint_c_from_T_RH_ice, rh_ice_percent, rh_percent,
)


def test_ice_monotonicity_and_inverse():
    temperatures = np.linspace(-100.0, 0.01, 10001)
    pressures = esat_ice_hpa(temperatures)
    assert np.all(np.diff(pressures) > 0)
    np.testing.assert_allclose(T_from_e_ice(pressures), temperatures, rtol=0, atol=1e-11)


def test_ice_triple_point_and_inverse_near_zero():
    np.testing.assert_allclose(esat_ice_hpa(0.01) * 100, 611.657, rtol=0, atol=1e-10)
    temperatures = np.array([-1e-6, -1e-9, 0.0, 1e-9, 1e-6, 0.01])
    np.testing.assert_allclose(T_from_e_ice(esat_ice_hpa(temperatures)), temperatures,
                               rtol=0, atol=1e-12)


@pytest.mark.parametrize('function', [esat_ice_hpa, dln_esat_ice_dT])
def test_ice_temperature_domain(function):
    temperatures = np.array([[-100, 0.01, -20],
        [np.nextafter(-100., -np.inf), np.nextafter(0.01, np.inf), np.nan],
        [np.inf, -np.inf, -coeffs_ice['b']]])
    with np.errstate(all='raise'):
        result = function(temperatures)
    assert result.shape == temperatures.shape
    assert np.isfinite(result[0]).all()
    assert np.isnan(result[1:]).all()


def test_ice_inverse_pressure_domain():
    lower, upper = esat_ice_hpa([-100., 0.01])
    pressures = np.array([[lower, upper, esat_ice_hpa(-20.)],
        [np.nextafter(lower, 0), np.nextafter(upper, np.inf), 0],
        [-1, np.inf, np.nan]])
    with np.errstate(all='raise'):
        result = T_from_e_ice(pressures)
    np.testing.assert_allclose(result[0], [-100, 0.01, -20], atol=1e-11)
    assert np.isnan(result[1:]).all()


def test_ice_derivative_matches_finite_difference():
    temperatures = np.linspace(-99, -0.1, 30)
    step = 1e-4
    numeric = (np.log(esat_ice_hpa(temperatures + step))
               - np.log(esat_ice_hpa(temperatures - step))) / (2*step)
    np.testing.assert_allclose(dln_esat_ice_dT(temperatures), numeric, rtol=1e-8)


def test_ice_relative_humidity_and_frostpoint():
    temperatures = np.array([[-80.], [-20.], [0.02]])
    rhs = np.array([0., 50., 100., 110., -1., np.inf, np.nan])
    with np.errstate(all='raise'):
        frost = frostpoint_c_from_T_RH_ice(temperatures, rhs)
    assert frost.shape == (3, 7)
    assert np.isnan(frost[:, 0]).all() and np.isnan(frost[:, 4:]).all()
    assert np.isnan(frost[2]).all()
    np.testing.assert_allclose(frost[:2, 2], temperatures[:2, 0], atol=1e-11)
    assert np.all(frost[:2, 1] < temperatures[:2, 0])
    assert np.all(frost[:2, 3] > temperatures[:2, 0])
    assert np.isnan(frostpoint_c_from_T_RH_ice(-100, 50))
    pressure = esat_ice_hpa(-20) * 1.1
    np.testing.assert_allclose(rh_ice_percent(-20, pressure), 110)
    assert rh_percent(-20, pressure) < 110  # RH over ice and water differ.
    assert esat_ice_hpa(-20) < esat_water_hpa(-20)
    np.testing.assert_allclose(rh_ice_percent(-20, [0, pressure, -1, np.inf, np.nan]),
                              [0, 110, np.nan, np.nan, np.nan])


def test_ice_empty_arrays_and_immutable_coefficients():
    for function in (esat_ice_hpa, dln_esat_ice_dT, T_from_e_ice):
        assert function(np.empty((0, 2))).shape == (0, 2)
    assert rh_ice_percent([], 1).shape == (0,)
    assert frostpoint_c_from_T_RH_ice([], 50).shape == (0,)
    with pytest.raises(TypeError):
        coeffs_ice['a'] = 1
    with pytest.raises(TypeError):
        coeffs_ice['domain_c']['min'] = -200


@pytest.mark.parametrize('lower,rmse_pct,max_pct', [
    (-40., 0.003138, 0.004325), (-100., 0.048709, 0.100000),
])
def test_ice_accuracy_against_iapws(lower, rmse_pct, max_pct):
    iapws = pytest.importorskip('iapws._iapws')
    temperatures = np.append(np.arange(lower, 0.01, 0.05), 0.01)
    reference = np.array([iapws._Sublimation_Pressure(float(t + 273.15))*1e4
                          for t in temperatures])
    err = (esat_ice_hpa(temperatures)/reference - 1)*100
    assert np.sqrt(np.mean(err**2)) < rmse_pct * 1.01
    assert np.max(abs(err)) < max_pct * 1.01
