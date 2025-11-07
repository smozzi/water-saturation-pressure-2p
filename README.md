# water-sat-pressure-2p

**A compact and accurate formula for the saturated vapor pressure of liquid water (including the supercooled regime), with a closed-form inverse.**

Fitted to both **IAPWS-95** for $`T \ge 0.01\,^\circ\mathrm{C}`$ and to **Murphy–Koop (2005)** for $`T < 0.01\,^\circ\mathrm{C}`$.

This work is derived from [**HarmoClimate**](https://github.com/smozzi/harmo-climate), a location-tuned climate-baseline project.

## Formula

```math
\ln E_s(T) = E_0 + \frac{a\,T}{b + T} + \frac{c\,T}{d + T} \qquad \text{with } T\text{ in }^\circ\mathrm{C},\ E_s\text{ in hPa.}
```

| Coefficient | Value |
| --- | --- |
| `E0` | 1.810270925564 |
| `a`  | 269.265582773152 |
| `b`  | 323.238664916362 |
| `c`  | −253.834491723435 |
| `d`  | 333.837330281331 |

*A closed-form inverse $`T(E_s)`$ is provided in the documentation.*

## Accuracy (relative error)

Percent-scale metrics from the latest benchmark:

**Primary range**
| Reference | Range (°C) | RMSE (\%) | Max abs. err (\%) |
| --- | --- | ---: | ---: |
| IAPWS-95 | 0.01 → 60 | **0.00012** | **0.00028** |
| Murphy–Koop | −25 → 0.01 | **0.00193** | **0.00339** |

**Extended range**
| Reference | Range (°C) | RMSE (\%) | Max abs. err (\%) |
| --- | --- | ---: | ---: |
| IAPWS-95 | 0.01 → 100 | 0.01101 | 0.04333 |
| Murphy–Koop | −40 → 0.01 | 0.08606 | 0.33960 |

*Scope:* liquid water only (including supercooled); no sublimation over ice coefficients provided for now.

## Installation
```bash
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -e .[dev]
```

## Quickstart
```python
import numpy as np
from wsp2p import esat_water_hpa, T_from_e_water, rh_percent, dewpoint_c_from_T_RH

T = np.linspace(-20.0, 40.0, 5)
Es = esat_water_hpa(T)             # saturated vapor pressure (hPa)
T_back = T_from_e_water(Es)        # exact inversion (°C)
RH = rh_percent(25.0, 18.0)        # ≈ 75 % for 25 °C & e = 18 hPa
dew = dewpoint_c_from_T_RH(30.0, 68.0)  # dew point (°C) without iteration
```

## API surface
| Function | Description |
| --- | --- |
| `esat_water_hpa(T_c)` | Saturation vapor pressure (hPa) over liquid + supercooled water using the two-pole log-pressure form. |
| `T_from_e_water(e_hpa)` | Closed-form inversion of the two-pole relation, returns °C and clamps to the domain −40…100 °C. |
| `rh_percent(T_c, e_hpa)` | Relative humidity (%) computed directly from saturation pressure; clipped to [0, 100]. |
| `dewpoint_c_from_T_RH(T_c, rh_percent)` | Algebraic dew point from ambient T/RH by chaining `esat` and the inverse. |
| `specific_humidity_kg_per_kg(T_c, rh_percent, p_hpa)` | Moist-air specific humidity using EPS = 0.621945 without iterative solvers. |
| `dln_esat_dT(T_c)` | Analytic derivative of ln Es for sensitivity work or adjoints. |
| `degc_to_kelvin`, `kelvin_to_degc`, `pa_to_hpa`, `hpa_to_pa` | Unit conversions (NumPy-broadcast friendly). |

## Domain & units
- Validated for `src/wsp2p/coeffs.json -> domain_c`: currently −40 to +100 °C.
- Inputs/outputs for vapor pressure use hPa; convert with helpers when necessary.
- Outputs are clipped to physically meaningful ranges (e ≥ 0, RH ∈ [0, 100], q ≥ 0).

## Benchmarks & figures
Run `notebooks/Benchmarks.ipynb` (NumPy, Matplotlib, iapws) to regenerate:
- `docs/figures/esat_curve.png` – Es(T) comparison vs references.
- `docs/figures/abs_error.png` – absolute error vs IAPWS-95 / Murphy–Koop.
- `docs/figures/rel_error.png` – relative error (%) with zoom near the triple point.
The notebook also prints a Markdown table with RMSE(log) and max |err| ready to paste into reports.

## Citations
- W. Wagner & A. Pruß (June 2002). *The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use*. Journal of Physical and Chemical Reference Data, **31**. doi:10.1063/1.1461829.
- D.M. Murphy & T. Koop (April 2005). *Review of the vapour pressures of supercooled water for atmospheric applications*. Quarterly Journal of the Royal Meteorological Society, **131**, 1539–1565. doi:10.1256/qj.04.94.

## License
- Code: Apache License (see `LICENSE`).
- Notebooks & rendered figures: CC-BY-4.0 attribution requested when reusing visuals.

## Repository layout
```
water-sat-pressure-2p/
├─ src/wsp2p/
│  ├─ coeffs.json
│  ├─ __init__.py
│  └─ esat.py
├─ tests/test_numeric.py
├─ notebooks/Benchmarks.ipynb
├─ docs/figures/
├─ CITATION.cff
├─ LICENSE
└─ pyproject.toml
```
