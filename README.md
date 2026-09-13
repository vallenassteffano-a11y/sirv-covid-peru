# SIRV model of COVID-19 vaccination impact in Peru

Extends the standard SIR epidemiological model with a vaccinated compartment (V)
to estimate how much Peru's vaccination programme reduced infections during the
third wave, 15 December 2021 – 16 March 2022.

Written for an IB Mathematics Extended Essay (2025).

## Results

Fitted parameters: β = 0.951, γ = 0.231

|                          | SIR (no vaccination) | SIRV (with vaccination) |
| ------------------------ | -------------------- | ----------------------- |
| Peak active cases        | 1,065,611             | 459,925                  |
| Total infections, 90 days| 8,302,462             | 4,460,496                |

Vaccination averted an estimated **3,841,966 infections (46.3% reduction)** and
lowered the peak by **605,686** cases.

![Model fit against real active-case data](fit_vs_data.png)
![SIR vs SIRV infected curves](sir_vs_sirv.png)
![Full SIRV compartment dynamics](sirv_simulation.png)

## Scripts

- `fit_parameters.py` — fits transmission (β) and recovery (γ) rates to the 37 real
  active-case observations using least-squares minimisation with SciPy's L-BFGS-B optimizer
- `sirv_simulation.py` — solves the four-compartment SIRV system via Euler's method (Δt = 0.1, 900 steps)
- `sir_vs_sirv.py` — runs both models on the same fitted parameters

Requires `numpy`, `scipy`, `matplotlib`.

## Data

- `active_cases_peru.csv` — 37 active-case observations (Worldometer)
- `vaccinations_peru.csv` — 91 days of daily second-dose administrations (Peruvian Ministry of Health)
- `vaccination_fit.png` — Excel trendline used to derive the degree-5 polynomial v(t)
  fed into the simulation as the vaccination rate

Initial conditions come from Peru's National Institute of Statistics (N = 33,350,300),
Worldometer case counts, and reported Pfizer/Sinopharm efficacy rates. t = 0 is set
at 15 December 2021, the start of the third wave, not the start of the pandemic —
so the model begins with a substantial pre-existing recovered and vaccinated
population rather than a fully susceptible one.

R0 (initial recovered) for the SIRV model is the cumulative case count from 15 days
before t0, reflecting an assumed 15-day recovery period. For the counterfactual SIR
model (no vaccination), R0 is that same figure *plus* V0 — everyone who had already
received a second dose before t0 is folded into "recovered" instead, since the
no-vaccination scenario needs somewhere to place people who are already immune for
either reason while keeping the total population N constant.

## Assumptions

**Inherited from the standard SIR model:**
- Total population N is constant — no births, deaths, or migration
- Contact between susceptible and infected individuals, and thus new infections, is
  proportional to S·I, scaled by transmission rate β and normalised by N
- Infected individuals recover at a constant rate γ (equivalent to a fixed average
  infectious period)
- No incubation period — individuals are infectious immediately upon infection

**Added for the SIRV extension:**
- Only second-dose recipients count as "vaccinated," to avoid double-counting people
  who received both doses
- Of second-dose recipients, only a fraction v_eff develop full immunity; that
  fraction becomes immune immediately (no delay for immune response to develop)
- Vaccinations are drawn only from the susceptible (S) compartment — the model
  assumes no already-infected or already-recovered individuals received a second
  dose during the simulation window
- Once immune via vaccination, individuals (like recovered individuals) cannot be
  reinfected or vaccinated again

## Limitations

Constant population, permanent immunity after infection or vaccination, immediate
immunity post-dose, and no latency period. Only 37 reported data points across 90
days, so the fitted curve underestimates the observed peak — visible in the first
figure above. Initial conditions (S0, R0, V0) rely on several simplifying
assumptions since the simulation starts mid-pandemic rather than at the initial
outbreak, which likely underestimates S0 by overestimating R0 and V0 — especially
given waning immunity, which the model does not account for.

## Notes

The optimizer setup was adapted from a Stack Overflow discussion on least-squares
SIR fitting; the Euler implementation follows a Plus Maths article.

The vaccination rate function v(t) is a degree-5 polynomial trendline fitted in
Excel to `vaccinations_peru.csv` (R² = 0.557, a moderate fit given daily reporting
noise); see `vaccination_fit.png`. Vaccine efficacy (v_eff = 0.95) is based on
Pfizer's efficacy alone for the ongoing simulation (t ≥ 0), even though Sinopharm
was also administered during this period — a simplification for tractability. A
separate blended Pfizer/Sinopharm efficacy (0.872) was used only once, to calculate
V0, the count of already-immune vaccinated individuals before t0.
