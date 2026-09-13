import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

beta = 0.951
gamma = 0.231
veff = 0.95
dt = 0.1
t_max = 90
num_steps = int(t_max / dt) + 1
t_values = np.linspace(0, t_max, num_steps)


def v(t):
    raw = -0.00115 * t**5 + 0.24004 * t**4 - 17.88265 * t**3 + 602.54558 * t**2 - 10052.19822 * t + 120899.04357
    return max(raw, 0)


S0_sir = 12933465
I0_sir = 19035
R0_sir = 20397800
N_sir = S0_sir + I0_sir + R0_sir

S_sir = np.zeros(num_steps)
I_sir = np.zeros(num_steps)
R_sir = np.zeros(num_steps)
S_sir[0], I_sir[0], R_sir[0] = S0_sir, I0_sir, R0_sir

for t in range(1, num_steps):
    dS = -beta * S_sir[t-1] * I_sir[t-1] / N_sir
    dI = beta * S_sir[t-1] * I_sir[t-1] / N_sir - gamma * I_sir[t-1]
    dR = gamma * I_sir[t-1]
    S_sir[t] = S_sir[t-1] + dS * dt
    I_sir[t] = I_sir[t-1] + dI * dt
    R_sir[t] = R_sir[t-1] + dR * dt

S0_sirv = 12933465
I0_sirv = 19035
R0_sirv = 2236351
V0_sirv = 18161449
N_sirv = S0_sirv + I0_sirv + R0_sirv + V0_sirv

S_sirv = np.zeros(num_steps)
I_sirv = np.zeros(num_steps)
R_sirv = np.zeros(num_steps)
V_sirv = np.zeros(num_steps)
S_sirv[0], I_sirv[0], R_sirv[0], V_sirv[0] = S0_sirv, I0_sirv, R0_sirv, V0_sirv

for t in range(1, num_steps):
    t_now = t_values[t-1]
    vt = veff * v(t_now)
    dS = -beta * S_sirv[t-1] * I_sirv[t-1] / N_sirv - vt
    dI = beta * S_sirv[t-1] * I_sirv[t-1] / N_sirv - gamma * I_sirv[t-1]
    dR = gamma * I_sirv[t-1]
    dV = vt
    S_sirv[t] = S_sirv[t-1] + dS * dt
    I_sirv[t] = I_sirv[t-1] + dI * dt
    R_sirv[t] = R_sirv[t-1] + dR * dt
    V_sirv[t] = V_sirv[t-1] + dV * dt

I_sir_m = I_sir / 1e6
I_sirv_m = I_sirv / 1e6

plt.figure(figsize=(10, 6))
plt.plot(t_values, I_sir_m, label="Infected (SIR)", color='red', linewidth=2)
plt.plot(t_values, I_sirv_m, label="Infected (SIRV)", color='purple', linewidth=2, linestyle='--')
plt.title("Comparison of Infected Population: SIR vs SIRV", fontsize=14)
plt.xlabel("Time(days)", fontsize=12)
plt.ylabel("Infected Individuals (Millions)", fontsize=12)
plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:.1f}M"))
plt.legend()
plt.grid(True)
plt.savefig('sir_vs_sirv.png', dpi=150, bbox_inches='tight')
