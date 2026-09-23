"""
Computes the headline results in the README and essay section 3.2:
peak active cases and total infections for the SIR (no new vaccination)
and SIRV (with vaccination) models, and the difference between them.
"""
import numpy as np

beta, gamma = 0.951, 0.231     # fitted in fit_parameters.py
veff = 0.95                    # Pfizer efficacy, used for t >= 0
dt, t_max = 0.1, 90
num_steps = int(t_max / dt) + 1
t_values = np.linspace(0, t_max, num_steps)


def v(t):
    """Daily second-dose rate: degree-5 polynomial fitted to vaccinations_peru.csv."""
    raw = (-0.00115 * t**5 + 0.24004 * t**4 - 17.88265 * t**3
           + 602.54558 * t**2 - 10052.19822 * t + 120899.04357)
    return max(raw, 0)


def simulate(S0, I0, R0, V0, vaccinate):
    """Euler's method for SIRV; with vaccinate=False this is the plain SIR model."""
    N = S0 + I0 + R0 + V0
    S, I, R, V = (np.zeros(num_steps) for _ in range(4))
    S[0], I[0], R[0], V[0] = S0, I0, R0, V0

    for n in range(1, num_steps):
        new_infections = beta * S[n-1] * I[n-1] / N
        vaccinated = veff * v(t_values[n-1]) if vaccinate else 0
        S[n] = S[n-1] + (-new_infections - vaccinated) * dt
        I[n] = I[n-1] + (new_infections - gamma * I[n-1]) * dt
        R[n] = R[n-1] + (gamma * I[n-1]) * dt
        V[n] = V[n-1] + vaccinated * dt

    return I, R


# Initial conditions (derived in initial_conditions.py)
I_sirv, R_sirv = simulate(12_933_465, 19_035, 2_236_351, 18_161_449, vaccinate=True)
I_sir, R_sir = simulate(12_933_465, 19_035, 20_397_800, 0, vaccinate=False)

# Total infections over the 90 days = I(90) + R(90) - R0  (essay section 3.2)
total_sirv = I_sirv[-1] + R_sirv[-1] - 2_236_351
total_sir = I_sir[-1] + R_sir[-1] - 20_397_800

averted = total_sir - total_sirv
peak_reduction = I_sir.max() - I_sirv.max()

print(f"{'':28}{'SIR (no new vaccination)':>26}{'SIRV (with vaccination)':>26}")
print(f"{'Peak active cases':28}{I_sir.max():>26,.0f}{I_sirv.max():>26,.0f}")
print(f"{'Total infections, 90 days':28}{total_sir:>26,.0f}{total_sirv:>26,.0f}")
print()
print(f"Infections averted: {averted:,.0f} ({averted / total_sir * 100:.1f}% reduction)")
print(f"Peak reduced by:    {peak_reduction:,.0f} active cases")
