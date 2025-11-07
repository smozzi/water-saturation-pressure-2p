"""Numerical regression tests for the two-pole saturation formulation."""

from __future__ import annotations

import math

import numpy as np
import pytest

from wsp2p import T_from_e_water, esat_water_hpa

GRID = np.arange(-40.0, 100.0001, 0.05, dtype=np.float64)

# Reference relative-error targets taken from docs/figures/benchmark_table.md (2025-11-14)
# (rows covering -40–0.01 °C and 0.01–100 °C) with a 5% safety margin.
BENCHMARK_MARGIN = 1.05
WARM_RMSE_PCT = 0.011010 * BENCHMARK_MARGIN
WARM_MAX_PCT = 0.043330 * BENCHMARK_MARGIN
SC_RMSE_PCT = 0.086061 * BENCHMARK_MARGIN
SC_MAX_PCT = 0.339598 * BENCHMARK_MARGIN
EPS = np.finfo(np.float64).eps


def test_inverse_accuracy():
    es = esat_water_hpa(GRID)
    back = T_from_e_water(es)
    delta = np.max(np.abs(back - GRID))
    assert delta < 1e-6, f"Inverse drift {delta:.3e} °C"


def test_monotonicity():
    es = esat_water_hpa(GRID)
    diffs = np.diff(es)
    assert np.all(diffs > 0.0), "Es(T) must be strictly increasing"


def test_continuity_near_freezing():
    for pivot in (-1e-4, 0.0, 0.01):
        offsets = np.array([-1e-6, 0.0, 1e-6], dtype=np.float64)
        sample = esat_water_hpa(pivot + offsets)
        span = float(np.ptp(sample))
        assert span < 1e-6, f"Continuity failed near T={pivot} °C"


def _iapws_available() -> bool:
    try:
        import iapws  # noqa: F401
    except Exception:
        return False
    return True


def _iapws_saturation_water(T_c: np.ndarray) -> np.ndarray:
    from iapws import IAPWS95

    T_k = T_c + 273.15
    T_k = T_k.clip(273.16, 373.15)
    out = np.full_like(T_c, np.nan, dtype=np.float64)
    vals = []
    for T_val in T_k:
        vals.append(IAPWS95(T=T_val, x=0).P * 1.0e4)  # MPa -> hPa
    out = np.array(vals, dtype=np.float64)
    return out


def _murphy_koop_water(T_c: np.ndarray) -> np.ndarray:
    T_k = T_c + 273.15
    ln_p = (
        54.842763
        - 6763.22 / T_k
        - 4.21 * np.log(T_k)
        + 0.000367 * T_k
        + np.tanh(0.0415 * (T_k - 218.8))
        * (
            53.878
            - 1331.22 / T_k
            - 9.44523 * np.log(T_k)
            + 0.014025 * T_k
        )
    )
    return np.exp(ln_p) / 100.0  # Pa -> hPa


@pytest.mark.skipif(not _iapws_available(), reason="iapws not installed")
def test_accuracy_against_references():
    warm_T = np.arange(0.01, 100.0001, 0.05, dtype=np.float64)
    warm_ref = _iapws_saturation_water(warm_T)
    warm_model = esat_water_hpa(warm_T)
    mask = np.isfinite(warm_ref)
    assert np.all(mask), "IAPWS reference returned NaNs in warm domain"
    warm_ref = warm_ref[mask]
    warm_model = warm_model[mask]
    warm_rel_err = (warm_model - warm_ref) / warm_ref * 100.0
    warm_rmse = math.sqrt(np.mean(warm_rel_err**2))
    warm_max = float(np.max(np.abs(warm_rel_err)))
    assert warm_rmse < WARM_RMSE_PCT
    assert warm_max < WARM_MAX_PCT

    cold_T = np.arange(-40.0, 0.0, 0.05, dtype=np.float64)
    cold_ref = _murphy_koop_water(cold_T)
    cold_model = esat_water_hpa(cold_T)
    cold_rel_err = (cold_model - cold_ref) / cold_ref * 100.0
    cold_rmse = math.sqrt(np.mean(cold_rel_err**2))
    cold_max = float(np.max(np.abs(cold_rel_err)))
    assert cold_rmse < SC_RMSE_PCT
    assert cold_max < SC_MAX_PCT
