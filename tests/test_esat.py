import numpy as np
import pytest

from wsp2p.esat import (
    coeffs,
    degc_to_kelvin,
    dewpoint_c_from_T_RH,
    dln_esat_dT,
    esat_water_hpa,
    hpa_to_pa,
    kelvin_to_degc,
    pa_to_hpa,
    rh_percent,
    specific_humidity_kg_per_kg,
    T_from_e_water,
)


def test_temperature_conversion_roundtrip():
    temps_c = np.array([-40.0, -10.0, 0.0, 37.5, 100.0])
    expected_k = temps_c + 273.15
    temps_k = degc_to_kelvin(temps_c)
    np.testing.assert_allclose(temps_k, expected_k, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(kelvin_to_degc(temps_k), temps_c, rtol=0.0, atol=1e-12)


def test_pressure_conversion_roundtrip():
    pressures_pa = np.array([50_000.0, 101_325.0, 150_000.0])
    hpa = pa_to_hpa(pressures_pa)
    np.testing.assert_allclose(hpa, pressures_pa / 100.0)
    np.testing.assert_allclose(hpa_to_pa(hpa), pressures_pa)


# Regression values independently evaluated with 50-digit Decimal arithmetic.
@pytest.mark.parametrize(
    "temp_c, expected",
    [
        (-40.0, 0.1897676661191643),
        (-20.0, 1.2550092325736795),
        (-5.0, 4.217693996937082),
        (0.0, 6.112114464992923),
        (15.0, 17.057945988932545),
        (30.0, 42.469717680585056),
        (60.0, 199.46402141608965),
        (100.0, 1013.7388627281676),
    ],
)

def test_esat_water_matches_regression_table(temp_c, expected):
    computed = esat_water_hpa(temp_c)
    np.testing.assert_allclose(computed, expected, rtol=1e-11)


def test_dln_esat_matches_finite_difference():
    temps = np.linspace(-35.0, 95.0, 25)
    analytic = dln_esat_dT(temps)
    eps = 1e-4
    ln_es_plus = np.log(esat_water_hpa(temps + eps))
    ln_es_minus = np.log(esat_water_hpa(temps - eps))
    numeric = (ln_es_plus - ln_es_minus) / (2.0 * eps)
    np.testing.assert_allclose(analytic, numeric, rtol=5e-5, atol=1e-6)


def test_T_from_e_water_known_pressures():
    e_values = np.array([0.5, 6.112103132923173, 50.0])
    expected_T = np.array([-30.200109657212924, -2.55130568176993e-05, 32.87426106343108])
    recovered = T_from_e_water(e_values)
    np.testing.assert_allclose(recovered, expected_T, rtol=0.0, atol=2e-6)


def test_T_from_e_water_invalid_inputs():
    invalid = np.array([0.0, -1.0, np.nan])
    out = T_from_e_water(invalid)
    assert np.isnan(out[0]) and np.isnan(out[1]) and np.isnan(out[2])


def test_inverse_is_stable_near_zero_temperature():
    temperatures = np.array([-1e-6, -1e-9, 0.0, 1e-9, 1e-6])
    recovered = T_from_e_water(esat_water_hpa(temperatures))
    np.testing.assert_allclose(recovered, temperatures, rtol=0.0, atol=1e-13)


def test_inverse_handles_zero_quadratic_coefficient(monkeypatch):
    import wsp2p.esat as module

    monkeypatch.setattr(module, 'coeffs', {**coeffs, 'c': 0.0})
    temperatures = np.linspace(-40.0, 100.0, 101)
    recovered = module.T_from_e_water(module.esat_water_hpa(temperatures))
    np.testing.assert_allclose(recovered, temperatures, rtol=0.0, atol=1e-12)


def test_rh_percent_behaves_expected():
    temp = np.array([22.0, 10.0, 5.0])
    e_inputs = np.array(
        [
            21.691850302240194,  # 82% relative humidity at 22 °C
            50.0,  # supersaturation must remain visible
            -5.0,  # invalid negative vapor pressure
        ]
    )
    expected = np.array([82.0, 50.0 / esat_water_hpa(10.0) * 100.0, np.nan])
    result = rh_percent(temp, e_inputs)
    np.testing.assert_allclose(result, expected)


def test_dewpoint_from_T_RH_regression():
    temp = np.array([30.0])
    rh = np.array([35.0])
    dew = dewpoint_c_from_T_RH(temp, rh)
    expected = np.array([12.880672917971367])
    np.testing.assert_allclose(dew, expected, atol=1e-6)


def test_specific_humidity_expected_value():
    temp = np.array([28.0])
    rh = np.array([65.0])
    pressure = np.array([950.0])
    expected = np.array([0.016257553502927877])
    q = specific_humidity_kg_per_kg(temp, rh, pressure)
    np.testing.assert_allclose(q, expected, rtol=0.0, atol=1e-12)


@pytest.mark.parametrize("function", [esat_water_hpa, dln_esat_dT])
def test_temperature_domain_is_enforced_elementwise(function):
    temps = np.array([[-40.0, 100.0, 0.0], [-40.01, 100.01, np.nan],
                      [np.inf, -np.inf, -coeffs["b"]]])
    with np.errstate(all="raise"):
        result = function(temps)
    assert result.shape == temps.shape
    assert np.isfinite(result[0]).all()
    assert np.isnan(result[1:]).all()


def test_inverse_rejects_pressures_outside_domain():
    e_min, e_max = esat_water_hpa([-40.0, 100.0])
    pressures = np.array([[e_min, e_max, esat_water_hpa(0.0)],
                          [np.nextafter(e_min, 0.0), np.nextafter(e_max, np.inf), 0.0],
                          [1e-300, 1e300, np.inf]])
    with np.errstate(all="raise"):
        result = T_from_e_water(pressures)
    np.testing.assert_allclose(result[0], [-40.0, 100.0, 0.0], atol=1e-12)
    assert np.isnan(result[1:]).all()


def test_dewpoint_does_not_clip_dry_air_or_invalid_rh():
    with np.errstate(all="raise"):
        result = dewpoint_c_from_T_RH(20.0, [0.0, 0.1, -1.0, np.nan, np.inf, 100.0, 110.0])
    assert np.isnan(result[:5]).all()
    np.testing.assert_allclose(result[5], 20.0, atol=1e-12)
    assert result[6] > 20.0


def test_specific_humidity_rejects_invalid_pressures():
    e = float(esat_water_hpa(20.0) * 0.5)
    pressures = [-10.0, 0.0, e, e / 2.0, np.inf, np.nan, 1000.0]
    with np.errstate(all="raise"):
        result = specific_humidity_kg_per_kg(20.0, 50.0, pressures)
    assert np.isnan(result[:6]).all()
    assert 0.0 < result[6] < 1.0
    assert np.isnan(specific_humidity_kg_per_kg(100.0, 100.0, 1000.0))


def test_specific_humidity_rh_validation_and_dry_air():
    with np.errstate(all="raise"):
        result = specific_humidity_kg_per_kg(20.0, [-1.0, np.nan, np.inf, 0.0, 110.0], 1000.0)
    assert np.isnan(result[:3]).all()
    assert result[3] == 0.0
    assert 0.0 < result[4] < 1.0


def test_broadcasting_with_invalid_elements():
    temps = np.array([[0.0], [20.0], [101.0]])
    rhs = np.array([50.0, 100.0, -1.0])
    es = esat_water_hpa(temps)
    recovered = T_from_e_water(es)
    np.testing.assert_allclose(recovered, temps * np.array([[1.0], [1.0], [np.nan]]))
    rh = rh_percent(temps, np.array([0.0, 5.0, -1.0]))
    assert rh.shape == (3, 3)
    assert np.isnan(rh[:, 2]).all() and np.isnan(rh[2]).all()
    dew = dewpoint_c_from_T_RH(temps, rhs)
    assert dew.shape == (3, 3)
    np.testing.assert_allclose(dew[:2, 1], temps[:2, 0], atol=1e-12)
    assert np.isnan(dew[:, 2]).all() and np.isnan(dew[2]).all()
    q = specific_humidity_kg_per_kg(temps, rhs, np.array([1000.0, 900.0, 0.0]))
    assert q.shape == (3, 3)
    assert np.isnan(q[:, 2]).all() and np.isnan(q[2]).all()
    assert np.isfinite(q[:2, :2]).all()


def test_coefficients_are_read_only_including_domain():
    with pytest.raises(TypeError):
        coeffs["a"] = 1.0
    with pytest.raises(TypeError):
        coeffs["domain_c"]["min"] = -100.0


def test_empty_arrays_are_supported():
    for function in (esat_water_hpa, dln_esat_dT, T_from_e_water):
        assert function(np.empty((0, 2))).shape == (0, 2)
    assert rh_percent([], 1.0).shape == (0,)
    assert dewpoint_c_from_T_RH([], 50.0).shape == (0,)
    assert specific_humidity_kg_per_kg([], 50.0, 1000.0).shape == (0,)
