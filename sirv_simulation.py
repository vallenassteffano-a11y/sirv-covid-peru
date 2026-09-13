import as np
import matplotlib.pyplot as plt
def v(t):
return -0.00115 * t**5 + 0.24004 * t**4 - 17.88265 * t**3 + 602.54558 * t**2 -
10052.19822 * t + 120899.04357
beta = 0.951
gamma = 0.231
veff = 0.95
N = 33350300
S0 = 12933465
I0 = 19035
R0 = 2236351
V0 = 18161449
t_max = 90
dt = 0.1
num_steps = int(t_max / dt) + 1
S = np.zeros(num_steps)
I = np.zeros(num_steps)
R = np.zeros(num_steps)
V = np.zeros(num_steps)
t_values = np.zeros(num_steps)
S[0] = S0
I[0] = I0
R[0] = R0
V[0] = V0
for t in range(1, num_steps):
t_values[t] = t * dt
dS = -beta * S[t-1] * I[t-1] / N - veff * v(t_values[t-1])
dI = beta * S[t-1] * I[t-1] / N - gamma * I[t-1]
dR = gamma * I[t-1]
dV = veff * v(t_values[t-1])
S[t] = S[t-1] + dS * dt
I[t] = I[t-1] + dI * dt
R[t] = R[t-1] + dR * dt
V[t] = V[t-1] + dV * dt
S_million = S / 1e6
I_million = I / 1e6
R_million = R / 1e6
V_million = V / 1e6
plt.figure(figsize=(10, 6))
plt.plot(t_values, S_million, label="Susceptible (S)", color='blue', linewidth=2)
plt.plot(t_values, I_million, label="Infected (I)", color='red', linewidth=2)
plt.plot(t_values, R_million, label="Recovered (R)", color='green', linewidth=2)
plt.plot(t_values, V_million, label="Vaccinated (V)", color='purple', linewidth=2)
plt.title("SIRV Model Simulation from December 15, 2021 to March 16, 2022",
fontsize=14)
plt.xlabel("Time(Days)", fontsize=12)
plt.ylabel("Number of individuals (Millions)", fontsize=12)
plt.legend()
plt.grid(True)
I_max = np.max(I_million)
I_at_90 = I_million[900]
R_at_90 = R_million[900]
R_at_0 = R_million[0]
print(f"Maximum Infected (I_max): {I_max}")
print(f"Infected at Day 90 (I(900)): {I_at_90}")
print(f"Recovered at Day 90 (R(900)): {R_at_90}")
print(f"Recovered at Initial Time (R(0)): {R_at_0}")
plt.tick_params(axis='y', labelsize=10)
plt.show()
