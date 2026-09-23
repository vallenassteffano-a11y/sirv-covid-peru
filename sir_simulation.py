"""
SIR model (no new vaccination) for Peru's third wave, solved with Euler's
method (essay Figure 15, Appendix 6).
"""
import numpy as np
import matplotlib.pyplot as plt

beta, gamma = 0.951, 0.231
dt, t_max = 0.1, 90
num_steps = int(t_max / dt) + 1
t_values = np.linspace(0, t_max, num_steps)

S0, I0, R0 = 12_933_465, 19_035, 20_397_800
N = S0 + I0 + R0

S, I, R = (np.zeros(num_steps) for _ in range(3))
S[0], I[0], R[0] = S0, I0, R0

for n in range(1, num_steps):
    new_infections = beta * S[n-1] * I[n-1] / N
    S[n] = S[n-1] - new_infections * dt
    I[n] = I[n-1] + (new_infections - gamma * I[n-1]) * dt
    R[n] = R[n-1] + (gamma * I[n-1]) * dt

print(f"Maximum infected (I_max): {I.max():,.0f}")
print(f"Infected at day 90:       {I[-1]:,.0f}")
print(f"Recovered at day 90:      {R[-1]:,.0f}")

plt.figure(figsize=(10, 6))
plt.plot(t_values, S / 1e6, label="Susceptible (S)", color="blue", linewidth=2)
plt.plot(t_values, I / 1e6, label="Infected (I)", color="red", linewidth=2)
plt.plot(t_values, R / 1e6, label="Recovered (R)", color="green", linewidth=2)
plt.title("SIR Model Simulation from December 15, 2021 to March 16, 2022", fontsize=14)
plt.xlabel("Time (days)", fontsize=12)
plt.ylabel("Number of individuals (Millions)", fontsize=12)
plt.legend()
plt.grid(True)
plt.savefig("sir_simulation.png", dpi=150, bbox_inches="tight")
plt.show()
