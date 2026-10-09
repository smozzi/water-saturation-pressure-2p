# API reference

The NumPy API provides separate quadratic saturation formulations for liquid
water and ice Ih, with closed-form inverses. Phase selection is explicit.

## Functions

| Function | Description |
| --- | --- |
| `esat_water_hpa(T_c)` | Saturation vapor pressure (hPa) over liquid + supercooled water using the four-coefficient log-pressure form. |
| `T_from_e_water(e_hpa)` | Closed-form inversion returning °C; pressures outside the domain return NaN. |
| `esat_ice_hpa(T_c)` | Saturation pressure over ice Ih (hPa), on −100…0.01 °C. |
| `T_from_e_ice(e_hpa)` | Closed-form frost point (°C); pressures outside the ice domain return NaN. |
| `rh_percent(T_c, e_hpa)` | Relative humidity (%); values above 100 retain supersaturation. |
| `rh_ice_percent(T_c, e_hpa)` | Relative humidity over ice (%), retaining supersaturation. |
| `dewpoint_c_from_T_RH(T_c, rh_percent_values)` | Algebraic dew point from ambient T/RH by chaining `esat` and the inverse. |
| `frostpoint_c_from_T_RH_ice(T_c, rh_ice_percent_values)` | Algebraic frost point from ambient T and RH referenced to ice. |
| `specific_humidity_kg_per_kg(T_c, rh_percent_values, p_hpa)` | Moist-air specific humidity using EPS = 0.621945 without iterative solvers. |
| `dln_esat_dT(T_c)` | Analytic derivative of ln Es for sensitivity work or adjoints. |
| `dln_esat_ice_dT(T_c)` | Analytic derivative of ice ln Es (K⁻¹). |
| `degc_to_kelvin`, `kelvin_to_degc`, `pa_to_hpa`, `hpa_to_pa` | Unit conversions (NumPy-broadcast friendly). |

## Domain and input behavior

Water functions cover **−40…100 °C**; ice functions cover **−100…0.01 °C**.
Inverse functions accept only pressures between the saturation pressures
at their respective domain endpoints. Invalid inputs return `NaN`
elementwise; values are not silently clipped into the domain.

- Compatible array shapes follow NumPy broadcasting; incompatible shapes raise `ValueError`.
- Nonfinite or out-of-domain temperatures produce `NaN`.
- Inverse pressure intervals are approximately 0.18977…1013.74 hPa for water and 0.0000140626…6.11657 hPa for ice.
- Relative humidity inputs must be finite and nonnegative. RH above 100% preserves supersaturation.
- Vapor pressure supplied to humidity functions must be finite and nonnegative.
- At RH = 0, dew/frost point is `NaN`. Points outside the relevant domain are also `NaN`, even if ambient temperature is valid.
- `specific_humidity_kg_per_kg` takes RH referenced to liquid water, and returns kg water per kg **moist air**. It requires finite total pressure greater than zero and greater than vapor pressure. At valid temperature and total pressure, RH = 0 gives zero specific humidity.
- Unit conversion helpers perform arithmetic only and do not enforce model domains.

```python
esat_water_hpa([-41.0, 20.0])       # [nan, valid saturation pressure]
dewpoint_c_from_T_RH(20.0, 0.1)    # nan: dew point is below −40 °C
```

## Inverse and derivative

For $y=\ln(e)-E_0$, solve $c z^2+a z-y=0$ using the increasing branch:

$$
z=\frac{2y}{a+\sqrt{a^2+4cy}},\qquad T=\frac{b z}{1-z}.
$$

This expression avoids cancellation near $y=0$ and also works when $c=0$.
Pressure inputs are checked against the saturation pressures at the domain
endpoints before solving. The analytic derivative is

$$
\frac{d\ln E_s}{dT}=(a+2cz)\frac{b}{(b+T)^2}.
$$

## General family and cubic alternative

The quadratic belongs to the family

$$
z=\frac{T}{b+T},\qquad \ln E_s(T)=E_0+\sum_{k=1}^{n} a_k z^k.
$$

Degree $n$ uses $n+2$ coefficients: $E_0$, $b$ and
$a_1,\ldots,a_n$. The package implements degree 2, with $a_1=a$ and $a_2=c$.

An optional degree-3 coefficient set improves agreement with IAPWS-95 over
0.01–100 °C. The default quadratic remains more accurate over 0.01–60 °C.
The cubic coefficients are provided separately in
[`coeffs_degree3.json`](coeffs_degree3.json); the package API uses
the quadratic fits above.

| Coefficient | Value |
| --- | ---: |
| `E0` | 1.8102749197 |
| `b` | 147.45412559 |
| `a1` | 10.715664579 |
| `a2` | 4.1962110938 |
| `a3` | 1.4533223391 |

| Reference | Range (°C) | RMSE (%) | Max abs. err (%) |
| --- | --- | ---: | ---: |
| IAPWS-95 | 0.01 → 100 | **0.000264** | **0.000376** |

These coefficients are rounded to 11 significant digits and evaluated only
on 0.01–100 °C.

## Citations

- W. Wagner & A. Pruß (June 2002). *The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use*. Journal of Physical and Chemical Reference Data, **31**. doi:10.1063/1.1461829.
- D.M. Murphy & T. Koop (April 2005). *Review of the vapour pressures of supercooled water for atmospheric applications*. Quarterly Journal of the Royal Meteorological Society, **131**, 1539–1565. doi:10.1256/qj.04.94.
- IAPWS (2011). *Revised Release on the Pressure along the Melting and Sublimation Curves of Ordinary Water Substance*, R14-08(2011). [Official release](https://iapws.org/technical-guidance/release/MeltSub).

[Return to README](../README.md).
