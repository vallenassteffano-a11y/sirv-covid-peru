import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
x_data = np.array([0, 2, 5, 7, 10, 12, 15, 17, 20, 23, 25, 28, 30, 33, 35,
38, 40, 43, 45, 48, 51, 53, 56, 58, 61, 63, 65, 67, 70, 72, 75, 78, 80, 83,
85, 88, 90])
y_data = np.array([19035, 19237, 19572, 19826, 24434, 24566, 30371, 31511,
40359, 63743, 94451, 148095, 215958, 321970, 393449, 518402, 563683, 607612,
590362, 538999, 437095, 387229, 277236, 228406, 187547, 157616, 127499,
113350, 81181, 66724, 52412, 33516, 29653, 23929, 22283, 18121, 15969])
S0 = 12933465
I0 = 19035
R0 = 2236351
V0 = 18161449
N = S0 + I0 + R0 + V0
v_eff = 0.95
def v(t):
return -0.00115 * t**5 + 0.24004 * t**4 - 17.88265 * t**3 + 602.54558 *
t**2 - 10052.19822 * t + 120899.04357
def euler_method(beta, gamma, dt, T):
t = np.arange(0, T + dt, dt)
S = np.zeros(len(t))
I = np.zeros(len(t))
R = np.zeros(len(t))
V = np.zeros(len(t))
S[0] = S0
I[0] = I0
R[0] = R0
V[0] = V0
for n in range(1, len(t)):
dSdt = -beta * S[n-1] * I[n-1] / N - v_eff * v(t[n-1])
dIdt = beta * S[n-1] * I[n-1] / N - gamma * I[n-1]
dRdt = gamma * I[n-1]
dVdt = v_eff * v(t[n-1])
S[n] = S[n-1] + dSdt * dt
I[n] = I[n-1] + dIdt * dt
R[n] = R[n-1] + dRdt * dt
V[n] = V[n-1] + dVdt * dt
return t, I
def objective(params):
beta, gamma = params
t, I_model = euler_method(beta, gamma, dt=0.1, T=90)
I_model_interpolated = np.interp(x_data, t, I_model)
return np.sum((I_model_interpolated - y_data) ** 2)
initial_guess = [0.3, 0.1]
result = minimize(objective, initial_guess, bounds=[(0, 10), (0, 1)])
beta_opt, gamma_opt = result.x
print(f"Optimized beta: {beta_opt}")
print(f"Optimized gamma: {gamma_opt}")
print(f"Sum of Least Squares Error: {result.fun}")
t_opt, I_opt = euler_method(beta_opt, gamma_opt, dt=0.1, T=90)
plt.figure(figsize=(10,6))
plt.plot(x_data, y_data, 'bo', label='Real Data (Active Cases)', markersize=6,
alpha=0.7)
plt.plot(np.arange(0, 90 + 0.1, 0.1), I_opt, 'r-', label='I(t)', linewidth=2)
45
plt.xlabel('Days')
plt.ylabel('Active Cases')
plt.legend()
plt.title('SIRV Model I(t) fit using Euler Method Solution and Least squares
method (LSM)')
plt.grid(True)
plt.show()
