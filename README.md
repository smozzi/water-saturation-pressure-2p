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

For $`y=\ln(e)-E_0`$, multiplying by both denominators gives

```math
A T^2+B T+C=0,\qquad
A=y-a-c,\quad B=y(b+d)-ad-cb,\quad C=ybd.
```

The implementation selects the root in the validated domain using
$`q=-\tfrac12(B+\operatorname{sign}(B)\sqrt{B^2-4AC})`$ and $`T=C/q`$
to avoid cancellation. Pressure inputs are checked against the saturation
pressures at the domain endpoints before solving.

## Accuracy (relative error)

Percent-scale metrics from the latest benchmark:

**Primary range**
| Reference | Range (°C) | RMSE (\%) | Max abs. err (\%) |
| --- | --- | ---: | ---: |
| IAPWS-95 | 0.01 → 60 | **0.00012** | **0.00028** |
| Murphy–Koop | −25 → 0.01 | **0.00192** | **0.00339** |

**Extended range**
| Reference | Range (°C) | RMSE (\%) | Max abs. err (\%) |
| --- | --- | ---: | ---: |
| IAPWS-95 | 0.01 → 100 | 0.01105 | 0.04345 |
| Murphy–Koop | −40 → 0.01 | 0.08601 | 0.33960 |

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
RH = rh_percent(25.0, 18.0)        # ≈ 57 % for 25 °C & e = 18 hPa
dew = dewpoint_c_from_T_RH(30.0, 68.0)  # dew point (°C) without iteration
```

## API surface
| Function | Description |
| --- | --- |
| `esat_water_hpa(T_c)` | Saturation vapor pressure (hPa) over liquid + supercooled water using the two-pole log-pressure form. |
| `T_from_e_water(e_hpa)` | Closed-form inversion returning °C; pressures outside the domain return NaN. |
| `rh_percent(T_c, e_hpa)` | Relative humidity (%); values above 100 retain supersaturation. |
| `dewpoint_c_from_T_RH(T_c, rh_percent_values)` | Algebraic dew point from ambient T/RH by chaining `esat` and the inverse. |
| `specific_humidity_kg_per_kg(T_c, rh_percent_values, p_hpa)` | Moist-air specific humidity using EPS = 0.621945 without iterative solvers. |
| `dln_esat_dT(T_c)` | Analytic derivative of ln Es for sensitivity work or adjoints. |
| `degc_to_kelvin`, `kelvin_to_degc`, `pa_to_hpa`, `hpa_to_pa` | Unit conversions (NumPy-broadcast friendly). |

## Domain & units
- Validated for `src/wsp2p/coeffs.json -> domain_c`: currently −40 to +100 °C.
- Inputs/outputs for vapor pressure use hPa; convert with helpers when necessary.
- Calculation functions return `NaN` **elementwise** for nonfinite, physically invalid, or out-of-domain inputs. Compatible array shapes follow NumPy broadcasting; incompatible shapes raise `ValueError`.
- Temperature inputs must lie in the inclusive domain −40…100 °C. The inverse accepts only pressures between `esat_water_hpa(-40)` and `esat_water_hpa(100)` (approximately 0.18976…1013.74 hPa).
- Vapor pressure for `rh_percent` must be finite and nonnegative. RH inputs must be finite and nonnegative; RH > 100 is permitted to represent supersaturation.
- At RH = 0, specific humidity is zero, while dew point is `NaN`. Dew points outside the domain also return `NaN`, even if the ambient temperature is valid.
- Specific humidity requires finite total pressure `p > 0` and vapor pressure `e < p`. It returns kg water per kg **moist air**, not the dry-air humidity ratio.
- Unit conversion helpers perform arithmetic only; they do not impose the saturation model's domain.
- The exported `coeffs` mapping and its `domain_c` mapping are read-only.

These rules replace the previous silent clipping of inverse temperatures, RH,
and total pressure. Callers that require clipping should apply it explicitly
before calling the API; callers should handle `NaN` results.

```python
esat_water_hpa([-41.0, 20.0])       # [nan, valid saturation pressure]
dewpoint_c_from_T_RH(20.0, 0.1)    # nan: dew point is below −40 °C
specific_humidity_kg_per_kg(20.0, 50.0, -10.0)  # nan: invalid total pressure
```

## Benchmarks & figures
Run `notebooks/Benchmarks.ipynb` (NumPy, Matplotlib, iapws) to regenerate:
- `docs/figures/esat_curve.png` – Es(T) comparison vs references.
- `docs/figures/abs_error.png` – absolute error vs IAPWS-95 / Murphy–Koop.
- `docs/figures/rel_error.png` – relative error (%) with zoom near the triple point.
The notebook includes exact range endpoints and prints a Markdown table with absolute (hPa) and relative (%) RMSE and maximum errors. These errors measure agreement with the references, not physical measurement uncertainty.

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
