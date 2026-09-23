"""
Sensitivity analysis of the SIRV model (essay section 1.8, Appendix 1).

Uses a hypothetical population (N = 10,000,000) and fixed beta, gamma to see
how the curves change when varying:
  - the vaccination rate v(t), with veff = 1          (essay Figures 3-5)
  - the vaccine efficacy veff, with v(t) = 20,000     (essay Figures 6-8)
"""
import numpy as np
import matplotlib.pyplot as plt

# Hypothetical initial conditions and parameters (essay Tables 3 and 4)
N = 10_000_000
S0, I0, R0, V0 = 4_000_000, 20_000, 3_000_000, 2_980_000
beta, gamma = 0.5, 0.1
days, dt = 100, 0.1


def simulate_sirv(v_rate, veff):
    """Solve the SIRV model with a constant vaccination rate using Euler's method."""
    time = np.arange(0, days, dt)
    S, I, R, V = (np.zeros(len(time)) for _ in range(4))
    S[0], I[0], R[0], V[0] = S0, I0, R0, V0

    for n in range(1, len(time)):
        new_infections = beta * S[n-1] * I[n-1] / N
        vaccinated = veff * v_rate
        S[n] = S[n-1] + (-new_infections - vaccinated) * dt
        I[n] = I[n-1] + (new_infections - gamma * I[n-1]) * dt
        R[n] = R[n-1] + (gamma * I[n-1]) * dt
        V[n] = V[n-1] + vaccinated * dt

    return time, S, I, R, V


scenarios = [
    # (v(t), veff, panel title)
    (10_000, 1.0, "Low vaccination rate\nv(t) = 10,000, veff = 1"),
    (20_000, 1.0, "Medium vaccination rate\nv(t) = 20,000, veff = 1"),
    (30_000, 1.0, "High vaccination rate\nv(t) = 30,000, veff = 1"),
    (20_000, 0.4, "Low vaccine efficacy\nv(t) = 20,000, veff = 0.4"),
    (20_000, 0.7, "Medium vaccine efficacy\nv(t) = 20,000, veff = 0.7"),
    (20_000, 1.0, "Perfect vaccine efficacy\nv(t) = 20,000, veff = 1"),
]

fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharey=True)

print(f"{'v(t)':>8} {'veff':>5} {'Peak infected':>15} {'Minimum S':>12}")
for ax, (v_rate, veff, title) in zip(axes.flat, scenarios):
    time, S, I, R, V = simulate_sirv(v_rate, veff)
    print(f"{v_rate:>8,} {veff:>5} {I.max():>15,.0f} {S.min():>12,.0f}")

    ax.plot(time, S / 1e6, label="Susceptible", color="blue")
    ax.plot(time, I / 1e6, label="Infected", color="red")
    ax.plot(time, R / 1e6, label="Recovered", color="green")
    ax.plot(time, V / 1e6, label="Vaccinated", color="purple")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("Time (days)")
    ax.grid(True)

axes[0, 0].set_ylabel("Population (millions)")
axes[1, 0].set_ylabel("Population (millions)")
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=11)
fig.suptitle("SIRV sensitivity analysis: varying v(t) (top) and veff (bottom)", fontsize=14)
fig.tight_layout(rect=(0, 0.05, 1, 1))

print("\nNote: at v(t) = 30,000 the susceptible population goes negative, because")
print("the model keeps removing veff * v(t) people per day even when S reaches 0.")

plt.savefig("sensitivity_analysis.png", dpi=150, bbox_inches="tight")
plt.show()
