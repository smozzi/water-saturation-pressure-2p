"""Two-pole saturated vapor pressure formulation for water."""

from __future__ import annotations

import json
from importlib.resources import files
from types import MappingProxyType
from typing import Any, Mapping

import numpy as np
from numpy.typing import ArrayLike

EPS = 0.621945
HPA = 100.0


def _load_coeffs() -> Mapping[str, Any]:
    coeffs_path = files("wsp2p").joinpath("coeffs.json")
    with coeffs_path.open("r", encoding="utf-8") as fh:
        data: dict[str, Any] = json.load(fh)
    required_scalars = {"E0", "a", "b", "c", "d"}
    missing = required_scalars - data.keys()
    if missing:
        missing_str = ", ".join(sorted(missing))
        raise KeyError(f"coeffs.json missing required keys: {missing_str}")
    for key in required_scalars:
        if isinstance(data[key], bool) or not isinstance(data[key], (int, float)):
            raise TypeError(f"Coefficient '{key}' must be numeric.")
        data[key] = float(data[key])
        if not np.isfinite(data[key]):
            raise ValueError(f"Coefficient '{key}' must be finite.")
    domain = data.get("domain_c")
    if not isinstance(domain, dict) or "min" not in domain or "max" not in domain:
        raise KeyError("coeffs.json must declare domain_c with 'min' and 'max'.")
    data["domain_c"] = {"min": float(domain["min"]), "max": float(domain["max"])}
    lower, upper = data["domain_c"]["min"], data["domain_c"]["max"]
    if not np.isfinite([lower, upper]).all() or lower >= upper:
        raise ValueError("domain_c must have finite, increasing bounds.")
    for key in ("b", "d"):
        if lower <= -data[key] <= upper:
            raise ValueError(f"Pole -{key} must lie outside domain_c.")
    data["domain_c"] = MappingProxyType(data["domain_c"])
    return MappingProxyType(data)


coeffs = _load_coeffs()


def _as_float_array(value: ArrayLike) -> np.ndarray:
    return np.asarray(value, dtype=np.float64)


def _validated_temperature(T_c: ArrayLike) -> np.ndarray:
    T = _as_float_array(T_c)
    valid = (
        np.isfinite(T)
        & (T >= coeffs["domain_c"]["min"])
        & (T <= coeffs["domain_c"]["max"])
    )
    return np.where(valid, T, np.nan)


def _vapor_pressure_from_rh(T_c: ArrayLike, rh_percent_values: ArrayLike) -> np.ndarray:
    rh = _as_float_array(rh_percent_values)
    rh = np.where(np.isfinite(rh) & (rh >= 0.0), rh, np.nan)
    with np.errstate(over="ignore"):
        e = esat_water_hpa(T_c) * (rh / 100.0)
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
    T = _validated_temperature(T_c)
    denom_b = coeffs["b"] + T
    denom_d = coeffs["d"] + T
    ln_es = coeffs["E0"] + (coeffs["a"] * T) / denom_b + (coeffs["c"] * T) / denom_d
    return np.exp(ln_es)


def dln_esat_dT(T_c: ArrayLike) -> np.ndarray:
    """Derivative of ln(Es); invalid/out-of-domain temperatures produce NaN."""
    T = _validated_temperature(T_c)
    term_a = coeffs["a"] * coeffs["b"] / ((coeffs["b"] + T) ** 2)
    term_c = coeffs["c"] * coeffs["d"] / ((coeffs["d"] + T) ** 2)
    return term_a + term_c


def _solve_quadratic(y: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=np.float64)

    a = coeffs["a"]
    b = coeffs["b"]
    c = coeffs["c"]
    d = coeffs["d"]

    # Multiplying y = a*T/(b+T) + c*T/(d+T) by both denominators
    # gives A*T**2 + B*T + C = 0. The other root is outside domain_c.
    a_coef = a + c
    A = y - a_coef
    B = y * (b + d) - (a * d + c * b)
    C = y * b * d

    # y is restricted to the validated pressure interval by the caller.
    disc = B * B - 4.0 * A * C
    sqrt_disc = np.sqrt(disc)
    sign_B = np.where(B >= 0.0, 1.0, -1.0)

    # Avoid subtracting nearly equal numbers; recover the physical root
    # through the product of the roots (C/A), including T=0 when y=0.
    q = -0.5 * (B + sign_B * sqrt_disc)

    T = C / q
    # Only correct floating-point drift at the endpoints, after input validation.
    return np.clip(T, coeffs["domain_c"]["min"], coeffs["domain_c"]["max"])


def T_from_e_water(e_hpa: ArrayLike) -> np.ndarray:
    """
    Closed-form inverse returning °C from vapor pressure in hPa.

    Nonfinite pressures and pressures outside the saturation-pressure interval
    corresponding to domain_c produce NaN, elementwise.
    """
    e = _as_float_array(e_hpa)
    out = np.full_like(e, np.nan, dtype=np.float64)
    e_min = esat_water_hpa(coeffs["domain_c"]["min"])
    e_max = esat_water_hpa(coeffs["domain_c"]["max"])
    valid = np.isfinite(e) & (e >= e_min) & (e <= e_max)
    if not np.any(valid):
        return out
    y = np.log(e[valid]) - coeffs["E0"]
    T_sol = _solve_quadratic(y)
    out[valid] = T_sol
    return out


def rh_percent(T_c: ArrayLike, e_hpa: ArrayLike) -> np.ndarray:
    """Relative humidity (%); supersaturation is retained as RH > 100.

    Negative/nonfinite vapor pressures or invalid temperatures produce NaN.
    """
    e = _as_float_array(e_hpa)
    e = np.where(np.isfinite(e) & (e >= 0.0), e, np.nan)
    with np.errstate(over="ignore"):
        rh = (e / esat_water_hpa(T_c)) * 100.0
    return np.where(np.isfinite(rh), rh, np.nan)


def dewpoint_c_from_T_RH(T_c: ArrayLike, rh_percent_values: ArrayLike) -> np.ndarray:
    """Liquid-water dew point (°C), with NumPy broadcasting.

    RH must be finite and nonnegative; supersaturation is permitted. A dew
    point outside domain_c, including RH=0, produces NaN rather than clipping.
    """
    return T_from_e_water(_vapor_pressure_from_rh(T_c, rh_percent_values))


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
        _vapor_pressure_from_rh(T_c, rh_percent_values), _as_float_array(p_hpa)
    )
    valid = np.isfinite(e) & np.isfinite(p) & (p > 0.0) & (e < p)
    # Work with the partial-pressure fraction to avoid overflow for large p.
    fraction = np.full(e.shape, np.nan, dtype=np.float64)
    np.divide(e, p, out=fraction, where=valid)
    return EPS * fraction / (1.0 - (1.0 - EPS) * fraction)
