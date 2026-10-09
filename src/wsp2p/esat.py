"""Quadratic saturated vapor pressure formulations for liquid water and ice Ih."""

from __future__ import annotations

import json
from importlib.resources import files
from types import MappingProxyType
from typing import Any, Mapping

import numpy as np
from numpy.typing import ArrayLike

EPS = 0.621945
HPA = 100.0


def _load_coeffs(filename: str = "coeffs.json") -> Mapping[str, Any]:
    coeffs_path = files("wsp2p").joinpath(filename)
    with coeffs_path.open("r", encoding="utf-8") as fh:
        data: dict[str, Any] = json.load(fh)
    required_scalars = {"E0", "a", "b", "c"}
    missing = required_scalars - data.keys()
    if missing:
        missing_str = ", ".join(sorted(missing))
        raise KeyError(f"{filename} missing required keys: {missing_str}")
    for key in required_scalars:
        if isinstance(data[key], bool) or not isinstance(data[key], (int, float)):
            raise TypeError(f"Coefficient '{key}' must be numeric.")
        data[key] = float(data[key])
        if not np.isfinite(data[key]):
            raise ValueError(f"Coefficient '{key}' must be finite.")
    domain = data.get("domain_c")
    if not isinstance(domain, dict) or "min" not in domain or "max" not in domain:
        raise KeyError(f"{filename} must declare domain_c with 'min' and 'max'.")
    data["domain_c"] = {"min": float(domain["min"]), "max": float(domain["max"])}
    lower, upper = data["domain_c"]["min"], data["domain_c"]["max"]
    if not np.isfinite([lower, upper]).all() or lower >= upper:
        raise ValueError("domain_c must have finite, increasing bounds.")
    if data["a"] <= 0.0 or data["b"] <= 0.0:
        raise ValueError("Coefficients 'a' and 'b' must be positive.")
    if lower <= -data["b"] <= upper:
        raise ValueError("Pole -b must lie outside domain_c.")
    endpoints = np.array([lower, upper])
    z = endpoints / (data["b"] + endpoints)
    if np.any(data["a"] + 2.0 * data["c"] * z <= 0.0):
        raise ValueError("Log-pressure must be strictly increasing on domain_c.")
    data["domain_c"] = MappingProxyType(data["domain_c"])
    return MappingProxyType(data)


coeffs = _load_coeffs()
coeffs_ice = _load_coeffs("coeffs_ice.json")


def _as_float_array(value: ArrayLike) -> np.ndarray:
    return np.asarray(value, dtype=np.float64)


def _validated_temperature(T_c: ArrayLike, co: Mapping[str, Any]) -> np.ndarray:
    T = _as_float_array(T_c)
    valid = (
        np.isfinite(T)
        & (T >= co["domain_c"]["min"])
        & (T <= co["domain_c"]["max"])
    )
    return np.where(valid, T, np.nan)


def _vapor_pressure_from_rh(
    T_c: ArrayLike, rh_percent_values: ArrayLike, co: Mapping[str, Any]
) -> np.ndarray:
    rh = _as_float_array(rh_percent_values)
    rh = np.where(np.isfinite(rh) & (rh >= 0.0), rh, np.nan)
    with np.errstate(over="ignore"):
        e = _esat_hpa(T_c, co) * (rh / 100.0)
    return np.where(np.isfinite(e), e, np.nan)


def degc_to_kelvin(T_c: ArrayLike) -> np.ndarray:
    """Convert °C to K."""
    return _as_float_array(T_c) + 273.15


def kelvin_to_degc(T_k: ArrayLike) -> np.ndarray:
    """Convert K to °C."""
    return _as_float_array(T_k) - 273.15


def pa_to_hpa(p_pa: ArrayLike) -> np.ndarray:
    """Pascal to hectopascal."""
    return _as_float_array(p_pa) / HPA


def hpa_to_pa(p_hpa: ArrayLike) -> np.ndarray:
    """Hectopascal to Pascal."""
    return _as_float_array(p_hpa) * HPA


def esat_water_hpa(T_c: ArrayLike) -> np.ndarray:
    """
    Saturated vapor pressure over liquid/supercooled water (hPa).

    Parameters
    ----------
    T_c : ArrayLike
        Temperature in degrees Celsius. Values outside coeffs["domain_c"]
        or nonfinite values produce NaN, elementwise.
    """
    return _esat_hpa(T_c, coeffs)


def esat_ice_hpa(T_c: ArrayLike) -> np.ndarray:
    """Saturation pressure over ice Ih (hPa), on −100…0.01 °C.

    Nonfinite and out-of-domain temperatures produce NaN, elementwise.
    """
    return _esat_hpa(T_c, coeffs_ice)


def _esat_hpa(T_c: ArrayLike, co: Mapping[str, Any]) -> np.ndarray:
    T = _validated_temperature(T_c, co)
    z = T / (co["b"] + T)
    ln_es = co["E0"] + z * (co["a"] + co["c"] * z)
    return np.exp(ln_es)


def dln_esat_dT(T_c: ArrayLike) -> np.ndarray:
    """Liquid-water dln(Es)/dT (K⁻¹); invalid temperatures produce NaN."""
    return _dln_esat_dT(T_c, coeffs)


def dln_esat_ice_dT(T_c: ArrayLike) -> np.ndarray:
    """Ice-Ih dln(Es)/dT (K⁻¹); invalid temperatures produce NaN."""
    return _dln_esat_dT(T_c, coeffs_ice)


def _dln_esat_dT(T_c: ArrayLike, co: Mapping[str, Any]) -> np.ndarray:
    T = _validated_temperature(T_c, co)
    denominator = co["b"] + T
    z = T / denominator
    return (co["a"] + 2.0 * co["c"] * z) * co["b"] / denominator**2


def _solve_quadratic(y: np.ndarray, co: Mapping[str, Any]) -> np.ndarray:
    y = np.asarray(y, dtype=np.float64)

    a = co["a"]
    b = co["b"]
    c = co["c"]
    # Solve c*z**2 + a*z - y = 0, where z = T/(b+T).
    # Rationalizing the increasing branch avoids cancellation near y=0;
    # this expression also works for c=0 without dividing by c.
    z = 2.0 * y / (a + np.sqrt(a * a + 4.0 * c * y))
    T = b * z / (1.0 - z)
    # Only correct floating-point drift at the endpoints, after input validation.
    return np.clip(T, co["domain_c"]["min"], co["domain_c"]["max"])


def T_from_e_water(e_hpa: ArrayLike) -> np.ndarray:
    """
    Closed-form inverse returning °C from vapor pressure in hPa.

    Nonfinite pressures and pressures outside the saturation-pressure interval
    corresponding to domain_c produce NaN, elementwise.
    """
    return _T_from_e(e_hpa, coeffs)


def T_from_e_ice(e_hpa: ArrayLike) -> np.ndarray:
    """Closed-form frost point (°C) from pressure (hPa), on −100…0.01 °C.

    Nonfinite pressures and pressures outside the ice saturation interval
    produce NaN, elementwise.
    """
    return _T_from_e(e_hpa, coeffs_ice)


def _T_from_e(e_hpa: ArrayLike, co: Mapping[str, Any]) -> np.ndarray:
    e = _as_float_array(e_hpa)
    out = np.full_like(e, np.nan, dtype=np.float64)
    e_min = _esat_hpa(co["domain_c"]["min"], co)
    e_max = _esat_hpa(co["domain_c"]["max"], co)
    valid = np.isfinite(e) & (e >= e_min) & (e <= e_max)
    if not np.any(valid):
        return out
    y = np.log(e[valid]) - co["E0"]
    T_sol = _solve_quadratic(y, co)
    out[valid] = T_sol
    return out


def rh_percent(T_c: ArrayLike, e_hpa: ArrayLike) -> np.ndarray:
    """Relative humidity over liquid water (%); supersaturation is retained.

    Negative/nonfinite vapor pressures or invalid temperatures produce NaN.
    """
    return _rh_percent(T_c, e_hpa, coeffs)


def rh_ice_percent(T_c: ArrayLike, e_hpa: ArrayLike) -> np.ndarray:
    """Relative humidity over ice (%), with supersaturation retained.

    Negative/nonfinite pressures or invalid temperatures produce NaN.
    """
    return _rh_percent(T_c, e_hpa, coeffs_ice)


def _rh_percent(T_c: ArrayLike, e_hpa: ArrayLike, co: Mapping[str, Any]) -> np.ndarray:
    e = _as_float_array(e_hpa)
    e = np.where(np.isfinite(e) & (e >= 0.0), e, np.nan)
    with np.errstate(over="ignore"):
        rh = (e / _esat_hpa(T_c, co)) * 100.0
    return np.where(np.isfinite(rh), rh, np.nan)


def dewpoint_c_from_T_RH(T_c: ArrayLike, rh_percent_values: ArrayLike) -> np.ndarray:
    """Liquid-water dew point (°C), with NumPy broadcasting.

    RH must be finite and nonnegative; supersaturation is permitted. A dew
    point outside domain_c, including RH=0, produces NaN rather than clipping.
    """
    return T_from_e_water(_vapor_pressure_from_rh(T_c, rh_percent_values, coeffs))


def frostpoint_c_from_T_RH_ice(T_c: ArrayLike, rh_ice_percent_values: ArrayLike) -> np.ndarray:
    """Frost point (°C) from temperature and RH referenced to ice, with broadcasting.

    RH must be finite and nonnegative; supersaturation is permitted. Frost
    points outside −100…0.01 °C, including RH=0, produce NaN.
    """
    return T_from_e_ice(_vapor_pressure_from_rh(T_c, rh_ice_percent_values, coeffs_ice))


def specific_humidity_kg_per_kg(
    T_c: ArrayLike,
    rh_percent_values: ArrayLike,
    p_hpa: ArrayLike,
) -> np.ndarray:
    """Specific humidity in kg water/kg moist air, with NumPy broadcasting.

    RH must be finite and nonnegative. Total pressure must be finite, positive,
    and strictly greater than vapor pressure. Invalid elements produce NaN.
    This is specific humidity q, not the dry-air humidity ratio W.
    """
    e, p = np.broadcast_arrays(
        _vapor_pressure_from_rh(T_c, rh_percent_values, coeffs), _as_float_array(p_hpa)
    )
    valid = np.isfinite(e) & np.isfinite(p) & (p > 0.0) & (e < p)
    # Work with the partial-pressure fraction to avoid overflow for large p.
    fraction = np.full(e.shape, np.nan, dtype=np.float64)
    np.divide(e, p, out=fraction, where=valid)
    return EPS * fraction / (1.0 - (1.0 - EPS) * fraction)
