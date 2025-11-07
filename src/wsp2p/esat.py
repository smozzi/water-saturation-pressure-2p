"""Two-pole saturated vapor pressure formulation for water."""

from __future__ import annotations

import json
from importlib.resources import files
from typing import Any, Dict

import numpy as np
from numpy.typing import ArrayLike

EPS = 0.621945
HPA = 100.0


def _load_coeffs() -> Dict[str, Any]:
    coeffs_path = files("wsp2p").joinpath("coeffs.json")
    with coeffs_path.open("r", encoding="utf-8") as fh:
        data: Dict[str, Any] = json.load(fh)
    required_scalars = {"E0", "a", "b", "c", "d"}
    missing = required_scalars - data.keys()
    if missing:
        missing_str = ", ".join(sorted(missing))
        raise KeyError(f"coeffs.json missing required keys: {missing_str}")
    for key in required_scalars:
        if not isinstance(data[key], (int, float)):
            raise TypeError(f"Coefficient '{key}' must be numeric.")
        data[key] = float(data[key])
    domain = data.get("domain_c")
    if not isinstance(domain, dict) or "min" not in domain or "max" not in domain:
        raise KeyError("coeffs.json must declare domain_c with 'min' and 'max'.")
    data["domain_c"] = {"min": float(domain["min"]), "max": float(domain["max"])}
    return data


coeffs = _load_coeffs()


def _as_float_array(value: ArrayLike) -> np.ndarray:
    return np.asarray(value, dtype=np.float64)


def _clip_temperature(T: np.ndarray) -> np.ndarray:
    return np.clip(T, coeffs["domain_c"]["min"], coeffs["domain_c"]["max"])


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
        Temperature in degrees Celsius. Valid for coeffs["domain_c"].
    """
    T = _as_float_array(T_c)
    denom_b = coeffs["b"] + T
    denom_d = coeffs["d"] + T
    ln_es = coeffs["E0"] + (coeffs["a"] * T) / denom_b + (coeffs["c"] * T) / denom_d
    es = np.exp(ln_es)
    # enforce positivity
    return np.maximum(es, 0.0)


def dln_esat_dT(T_c: ArrayLike) -> np.ndarray:
    """Derivative of ln(Es) with respect to temperature."""
    T = _as_float_array(T_c)
    term_a = coeffs["a"] * coeffs["b"] / ((coeffs["b"] + T) ** 2)
    term_c = coeffs["c"] * coeffs["d"] / ((coeffs["d"] + T) ** 2)
    return term_a + term_c


def _solve_quadratic(y: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=np.float64)

    a = coeffs["a"]
    b = coeffs["b"]
    c = coeffs["c"]
    d = coeffs["d"]

    a_coef = a + c
    A = y - a_coef
    B = y * (b + d) - (a * d + c * b)
    C = y * b * d

    disc = np.maximum(B * B - 4.0 * A * C, 0.0)
    sqrt_disc = np.sqrt(disc)
    sign_B = np.where(B >= 0.0, 1.0, -1.0)

    q = -0.5 * (B + sign_B * sqrt_disc)

    T = C / q
    return _clip_temperature(T)


def T_from_e_water(e_hpa: ArrayLike) -> np.ndarray:
    """
    Closed-form inverse that returns temperature (°C) from vapor pressure (hPa).
    """
    e = _as_float_array(e_hpa)
    out = np.full_like(e, np.nan, dtype=np.float64)
    valid = np.isfinite(e) & (e > 0.0)
    if not np.any(valid):
        return out
    y = np.log(e[valid]) - coeffs["E0"]
    T_sol = _solve_quadratic(y)
    out[valid] = T_sol
    return out


def rh_percent(T_c: ArrayLike, e_hpa: ArrayLike) -> np.ndarray:
    """Relative humidity (%) computed without iteration."""
    T = _as_float_array(T_c)
    e = _as_float_array(e_hpa)
    es = esat_water_hpa(T)
    rh = np.where(es > 0.0, (e / es) * 100.0, 0.0)
    return np.clip(rh, 0.0, 100.0)


def dewpoint_c_from_T_RH(T_c: ArrayLike, rh_percent_values: ArrayLike) -> np.ndarray:
    """Return dewpoint °C from ambient temperature and RH%."""
    T = _as_float_array(T_c)
    rh = np.clip(_as_float_array(rh_percent_values), 0.0, 100.0)
    es = esat_water_hpa(T)
    e = es * (rh / 100.0)
    return T_from_e_water(e)


def specific_humidity_kg_per_kg(
    T_c: ArrayLike,
    rh_percent_values: ArrayLike,
    p_hpa: ArrayLike,
) -> np.ndarray:
    """
    Specific humidity (kg/kg of moist air) using EPS ratio, no iteration required.
    """
    T = _as_float_array(T_c)
    rh = np.clip(_as_float_array(rh_percent_values), 0.0, 100.0)
    p = np.maximum(_as_float_array(p_hpa), 1.0)  # avoid zero/negative pressure
    e = esat_water_hpa(T) * (rh / 100.0)
    denom = p - (1.0 - EPS) * e
    q = EPS * e / denom
    return np.clip(q, 0.0, None)
