"""
Derives the initial conditions for the SIRV and SIR models at t = 0
(15 December 2021), following sections 2.3-2.5 of the essay.
"""

# Total population of Peru, 2021 (INEI)
N = 33_350_300

# Active cases on 15 December 2021 (Worldometer)
I0 = 19_035

# Second doses administered before t = 0
second_doses_before_t0 = 20_827_350

# Vaccine efficacy before t = 0: Pfizer and Sinopharm weighted equally,
# since the exact split of doses is unknown
pfizer_eff = 0.95
sinopharm_eff = 0.7934
veff_before_t0 = round((pfizer_eff + sinopharm_eff) / 2, 3)  # 0.872

# Only the fraction veff_before_t0 of second-dose recipients counts as immune
V0 = round(veff_before_t0 * second_doses_before_t0)

# Recovered: cumulative cases on 30 November 2021, assuming a 15-day
# recovery period (closest available date to 15 days before t = 0)
R0 = 2_236_351

# Susceptible: everyone else
S0 = N - I0 - R0 - V0

# SIR (no-vaccination) scenario: the already-vaccinated are folded into R
# so that N stays constant
R0_sir = R0 + V0

print(f"Blended vaccine efficacy before t0: {veff_before_t0}")
print()
print("SIRV initial conditions")
print(f"  S0 = {S0:,}")
print(f"  I0 = {I0:,}")
print(f"  R0 = {R0:,}")
print(f"  V0 = {V0:,}")
print(f"  N  = {S0 + I0 + R0 + V0:,}")
print()
print("SIR initial conditions (no vaccination)")
print(f"  S0 = {S0:,}")
print(f"  I0 = {I0:,}")
print(f"  R0 = {R0_sir:,}")
print(f"  N  = {S0 + I0 + R0_sir:,}")

# Check these match the values hard-coded in the other scripts
assert (S0, I0, R0, V0) == (12_933_465, 19_035, 2_236_351, 18_161_449)
assert R0_sir == 20_397_800
