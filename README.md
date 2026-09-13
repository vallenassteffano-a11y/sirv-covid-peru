# SIRV model of COVID-19 vaccination impact in Peru

Extends the standard SIR epidemiological model with a vaccinated compartment (V)
to estimate how much Peru's vaccination programme reduced infections during the
third wave, 15 December 2021 – 16 March 2022.

Written for an IB Mathematics Extended Essay (2025).

## Results

Fitted parameters: β = 0.951, γ = 0.231

|                          | SIR (no vaccination) | SIRV (with vaccination) |
| ------------------------ | -------------------- | ----------------------- |
| Peak active cases        | 1,065,611            | 459,925                 |
| Total infections, 90 days| 8,302,462            | 4,460,496               |

Vaccination averted an estimated **3,841,966 infections (46.3% reduction)** and
lowered the peak by **605,686** cases.

![Model fit against real active-case data](fit_vs_data.png)
![SIR vs SIRV infected curves](sir_vs_sirv.png)
![Full SIRV compartment dynamics](sirv_simulation.png)

## Scripts

- `fit_parameters.py` — fits transmission (β) and recovery (γ) rates to the 37 real
  active-case observations using least-squares minimisation with SciPy's L-BFGS-B optimizer
- `sirv_simulation.py` — solves the four-compartment SIRV system via Euler's method (Δt = 0.1, 900 steps)
- `sir_vs_sirv.py` — runs both models on the same fitted parameters and reports the difference

Requires `numpy`, `scipy`, `matplotlib`.

## Data

- `active_cases_peru.csv` — 37 active-case observations (Worldometer)
- `vaccinations_peru.csv` — 91 days of daily second-dose administrations (Peruvian Ministry of Health)

Initial conditions come from Peru's National Institute of Statistics (N = 33,350,300),
Worldometer case counts, and reported Pfizer/Sinopharm efficacy rates.

## Limitations

Constant population, permanent immunity after infection or vaccination, immediate
immunity post-dose, and no latency period. Only 37 reported data points across 90
days, so the fitted curve underestimates the observed peak — visible in the first
figure above.

## Notes

The optimizer setup was adapted from a Stack Overflow discussion on least-squares
SIR fitting; the Euler implementation follows a Plus Maths article.
