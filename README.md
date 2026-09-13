# SIRV model of COVID-19 vaccination impact in Peru

Extends the standard SIR epidemiological model with a vaccinated compartment (V)
to estimate how much Peru's vaccination programme reduced infections during the
third wave, 15 December 2021 – 16 March 2022.

## Files

- `fit_parameters.py` — fits transmission (β) and recovery (γ) rates to 37 real
  active-case observations using least-squares minimisation with SciPy's L-BFGS-B optimizer
- `sirv_simulation.py` — solves the SIRV system numerically via Euler's method (Δt = 0.1, 900 steps)
- `active_cases_peru.csv` — 37 active-case observations (Worldometer)
- `vaccinations_peru.csv` — 91 days of daily second-dose administrations (Peruvian Ministry of Health)

## Results

Fitted parameters: β = 0.951, γ = 0.231

Vaccination averted an estimated 3.84M infections (46.3% reduction) and lowered
peak active cases by 605,686.

![Model fit](fit_vs_data.png)
![SIR vs SIRV infected curves](sir_vs_sirv.png)

## Limitations

Constant population, permanent immunity after infection or vaccination, immediate
immunity post-dose, and no latency period. Only 37 reported data points across 90
days, so the fitted curve underestimates the observed peak.

## Notes

Written for an IB Mathematics Extended Essay (2025). The optimizer setup was adapted
from a Stack Overflow discussion on least-squares SIR fitting; the Euler implementation
follows a Plus Maths article.
