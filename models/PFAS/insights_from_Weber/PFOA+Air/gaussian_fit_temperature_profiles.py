# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Global Fit of Temperature Profiles — ODR
#
# Uses orthogonal distance regression (scipy.odr) to account for
# uncertainty in both x (distance, ~2 cm) and y (temperature, ~50 C).
#
# 13 shared parameters across all 7 profiles.
#
# Model:
#
# T(x, nominal) = baseline
#     + a1*(nominal - b1) * [G(x, 30-mu1, sigma1) + G(x, 30+mu1, sigma1)]
#     + a2*(nominal - b2) * [G(x, 30-mu2, sigma2) + G(x, 30+mu2, sigma2)]
#     + a3*(nominal - b3) * G(x, 30, sigma3)
#     + a4*(nominal - b4) * G(x, mu4, sigma4)

# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.odr import ODR, Model, RealData
from scipy.optimize import least_squares

# %%
# Load data
temperature_profiles = pd.read_csv('temperature_profiles_from_literature.csv')

nominal_temperatures = list(range(800, 1150, 50))

data = {}
for temp in nominal_temperatures:
    index_of_x = temperature_profiles.columns.get_loc(str(temp))
    index_of_y = index_of_x + 1

    distances_str = temperature_profiles.iloc[:, index_of_x]
    distances = np.array([float(x) for x in distances_str[1:]])

    temperatures_str = temperature_profiles.iloc[:, index_of_y]
    temperatures_at_distances = np.array([float(x) for x in temperatures_str[1:]])

    order = np.argsort(distances)
    data[temp] = (distances[order], temperatures_at_distances[order])

# %% [markdown]
# ## Model definition

# %%
MIDPOINT = 30.0

def sum_of_gaussians_odr(params, x_and_nominal):
    """ODR model function.

    scipy.odr passes the independent variables as a 2-row array:
      x_and_nominal[0] = x (distance)
      x_and_nominal[1] = nominal temperature
    """
    (baseline,
     a1, b1, mu1, sigma1,
     a2, b2, mu2, sigma2,
     a3, b3, sigma3,
     a4, b4, mu4, sigma4) = params

    x = x_and_nominal[0]
    nominal = x_and_nominal[1]
    mid = MIDPOINT

    return (baseline
            + a1 * (nominal - b1) * np.exp(-0.5 * ((x - (mid - mu1)) / sigma1) ** 2)
            + a1 * (nominal - b1) * np.exp(-0.5 * ((x - (mid + mu1)) / sigma1) ** 2)
            + a2 * (nominal - b2) * np.exp(-0.5 * ((x - (mid - mu2)) / sigma2) ** 2)
            + a2 * (nominal - b2) * np.exp(-0.5 * ((x - (mid + mu2)) / sigma2) ** 2)
            + a3 * (nominal - b3) * np.exp(-0.5 * ((x - mid) / sigma3) ** 2)
            + a4 * (nominal - b4) * np.exp(-0.5 * ((x - mu4) / sigma4) ** 2))

# %% [markdown]
# ## Get good initial parameters via ordinary least squares, then refine with ODR

# %%
# Stack all data
x_all = []
y_all = []
nominal_all = []
for temp in nominal_temperatures:
    xd, yd = data[temp]
    x_all.append(xd)
    y_all.append(yd)
    nominal_all.append(np.full_like(xd, temp))

x_all = np.concatenate(x_all)
y_all = np.concatenate(y_all)
nominal_all = np.concatenate(nominal_all)

# --- First: ordinary least squares to get a good starting point ---
def residuals_ls(params):
    (baseline,
     a1, b1, mu1, sigma1,
     a2, b2, mu2, sigma2,
     a3, b3, sigma3,
     a4, b4, mu4, sigma4) = params
    xn = np.vstack([x_all, nominal_all])
    y_pred = sum_of_gaussians_odr(params, xn)
    return y_pred - y_all

lower = np.array([0,    0.001, -5000,  5, 1.5,
                        0.001, -5000,  0, 1.5,
                        0.001, -5000,     1.5,
                        0.001, -5000,  8, 1.5])
upper = np.array([80,   20,     800, 28, 15,
                        20,     800, 14, 15,
                        20,     800,     15,
                        20,     800, 52, 15])

best_cost = np.inf
best_p0 = None

initial_guesses = [
    [25,  1.0, 700, 10, 5,  0.3, -1000, 4, 4,  0.5, -500, 3,  0.3,  500, 15,  3],
    [25,  0.5, 400, 12, 4,  0.2, -2000, 3, 5,  0.8, -200, 4,  0.5,  300, 14,  4],
    [25,  2.0, 750, 10, 6,  0.1, -3000, 4, 4,  0.3,-1000, 5,  0.2,  600, 16,  3],
    [30,  0.8, 500, 15, 5,  0.5,  -500, 5, 4,  1.0,    0, 4,  0.1,  200, 13,  5],
    [20,  1.5, 780,  9, 5,  0.05,-4000, 4, 4,  0.2,-2000, 4,  0.5,  700, 15,  2],
    [25,  0.3, 200, 12, 4,  0.4,  -800, 3, 5,  0.6, -300, 3,  0.4,  400, 16,  4],
    [25,  1.2, 600, 11, 5,  0.15,-2500, 4, 4,  0.4,-1500, 5,  0.3,  500, 14,  3],
    [25,  0.7, 300, 13, 4,  0.3, -1500, 3, 4,  0.9, -100, 3,  0.2,  350, 17,  3],
    [25,  1.0, 700, 10, 5,  0.1, -2000, 4, 4,  0.5, -200, 4,  0.5,  700, 42,  4],
    [25,  1.5, 780, 10, 5,  0.05,-3000, 4, 4,  0.4, -500, 4,  1.0,  750, 45,  5],
    [25,  0.8, 600, 11, 5,  0.2, -1000, 5, 3,  0.6, -100, 3,  0.8,  650, 40,  5],
    [25,  1.0, 700, 10, 5,  0.1, -2000, 4, 4,  0.5, -200, 4,  0.3,  600, 35,  8],
    [25,  1.2, 750, 10, 5,  0.08,-2500, 4, 4,  0.4, -300, 4,  0.5,  700, 38,  6],
]

print("Finding good starting parameters via least squares...")
for i, p0 in enumerate(initial_guesses):
    try:
        res = least_squares(residuals_ls, p0, bounds=(lower, upper),
                            method='trf', max_nfev=500000)
        if res.cost < best_cost:
            best_cost = res.cost
            best_p0 = res.x.copy()
            print(f"  Guess {i}: cost={res.cost:.1f} ** new best **")
        else:
            print(f"  Guess {i}: cost={res.cost:.1f}")
    except Exception as e:
        print(f"  Guess {i}: failed ({e})")

print(f"\nBest LS cost = {best_cost:.1f}")
print(f"Starting ODR from these parameters.\n")

# %% [markdown]
# ## Run ODR
#
# Uncertainties: sx = 2 cm in distance, sy = 50 C in temperature.
# The nominal temperature coordinate has no error (it's a known setpoint),
# so we set its uncertainty to a tiny value.

# %%
sx_dist = 2.0    # cm uncertainty in distance
sy_temp = 50.0   # C uncertainty in temperature
sx_nom = 1e-10   # nominal temperature is known exactly

# Build the 2-row independent variable array and matching uncertainties
xn_all = np.vstack([x_all, nominal_all])
sx_all = np.vstack([np.full_like(x_all, sx_dist),
                     np.full_like(nominal_all, sx_nom)])

odr_data = RealData(xn_all, y_all, sx=sx_all, sy=np.full_like(y_all, sy_temp))
odr_model = Model(sum_of_gaussians_odr)

odr_instance = ODR(odr_data, odr_model, beta0=best_p0, maxit=5000)
odr_result = odr_instance.run()

print("ODR result:")
odr_result.pprint()

(baseline,
 a1, b1, mu1, sigma1,
 a2, b2, mu2, sigma2,
 a3, b3, sigma3,
 a4, b4, mu4, sigma4) = odr_result.beta

print(f"\nFitted parameters:")
print(f"  baseline = {baseline:.2f} C")
print(f"  Outer pair:     a1={a1:.4f}, b1={b1:.1f}, mu1={mu1:.2f} cm, sigma1={sigma1:.2f} cm")
print(f"                  centers at {MIDPOINT-mu1:.1f} and {MIDPOINT+mu1:.1f} cm")
print(f"  Inner pair:     a2={a2:.4f}, b2={b2:.1f}, mu2={mu2:.2f} cm, sigma2={sigma2:.2f} cm")
print(f"                  centers at {MIDPOINT-mu2:.1f} and {MIDPOINT+mu2:.1f} cm")
print(f"  Center fill:    a3={a3:.4f}, b3={b3:.1f}, sigma3={sigma3:.2f} cm  (at {MIDPOINT} cm)")
print(f"  Asym. shoulder: a4={a4:.4f}, b4={b4:.1f}, mu4={mu4:.2f} cm, sigma4={sigma4:.2f} cm")

# Per-profile RMSE (in T only, for comparison with earlier fits)
print()
for temp in nominal_temperatures:
    xd, yd = data[temp]
    xn = np.vstack([xd, np.full_like(xd, temp)])
    y_pred = sum_of_gaussians_odr(odr_result.beta, xn)
    rmse = np.sqrt(np.mean((y_pred - yd) ** 2))
    print(f"  T={temp} C:  RMSE = {rmse:.1f} C")

# %% [markdown]
# ## Plot: all profiles overlay

# %%
colors = {800: '#4472C4',
          850: '#EE7D31',
          900: '#A5A5A5',
          950: '#FFC002',
          1000: '#5B9CD5',
          1050: '#70AE47',
          1100: '#264479'}

x_fine = np.linspace(0, 60, 300)

fig, ax = plt.subplots(figsize=(10, 6))

for temp in nominal_temperatures:
    xd, yd = data[temp]
    color = colors[temp]

    # data with error bars
    ax.errorbar(xd, yd, xerr=sx_dist, yerr=sy_temp,
                fmt='o', color=color, markersize=3, alpha=0.6,
                elinewidth=0.5, capsize=0)

    # fit line
    xn_fine = np.vstack([x_fine, np.full_like(x_fine, temp)])
    y_fit = sum_of_gaussians_odr(odr_result.beta, xn_fine)
    ax.plot(x_fine, y_fit, '-', color=color, linewidth=1.5, label=f'{temp} C')

ax.set_xlabel('Distance into reactor (cm)', fontweight='bold')
ax.set_ylabel('Temperature (C)', fontweight='bold')
ax.set_title('Global fit (ODR): 13 parameters, $\\sigma_x$=2 cm, $\\sigma_T$=50 C',
             fontweight='bold')
ax.legend(loc='upper left')
ax.set_xlim(0, 60)
ax.set_ylim(0, None)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontweight('bold')

plt.tight_layout()
plt.savefig('global_fit_odr_all.png', dpi=150, bbox_inches='tight')
plt.show()

# %% [markdown]
# ## Individual profile panels with component Gaussians

# %%
fig, axes = plt.subplots(3, 3, figsize=(15, 12), sharex=True)
axes_flat = axes.flatten()

for idx, temp in enumerate(nominal_temperatures):
    ax = axes_flat[idx]
    xd, yd = data[temp]
    color = colors[temp]

    ax.errorbar(xd, yd, xerr=sx_dist, yerr=sy_temp,
                fmt='o', color='black', markersize=3, alpha=0.6,
                elinewidth=0.5, capsize=0, label='data')

    xn_fine = np.vstack([x_fine, np.full_like(x_fine, temp)])
    y_fit = sum_of_gaussians_odr(odr_result.beta, xn_fine)
    ax.plot(x_fine, y_fit, '-', color=color, linewidth=2, label='ODR fit')

    # Component Gaussians (shifted up by baseline for display)
    b_off = baseline
    g_outer_l  = b_off + a1*(temp-b1)*np.exp(-0.5*((x_fine-(MIDPOINT-mu1))/sigma1)**2)
    g_outer_r  = b_off + a1*(temp-b1)*np.exp(-0.5*((x_fine-(MIDPOINT+mu1))/sigma1)**2)
    g_inner_l  = b_off + a2*(temp-b2)*np.exp(-0.5*((x_fine-(MIDPOINT-mu2))/sigma2)**2)
    g_inner_r  = b_off + a2*(temp-b2)*np.exp(-0.5*((x_fine-(MIDPOINT+mu2))/sigma2)**2)
    g_center   = b_off + a3*(temp-b3)*np.exp(-0.5*((x_fine-MIDPOINT)/sigma3)**2)
    g_shoulder = b_off + a4*(temp-b4)*np.exp(-0.5*((x_fine-mu4)/sigma4)**2)

    ax.plot(x_fine, g_outer_l, '--', color='gray', lw=1, alpha=0.5,
            label=f'outer ({MIDPOINT-mu1:.0f} cm)')
    ax.plot(x_fine, g_outer_r, '--', color='gray', lw=1, alpha=0.5)
    ax.plot(x_fine, g_inner_l, '-.', color='gray', lw=1, alpha=0.5,
            label=f'inner ({MIDPOINT-mu2:.0f} cm)')
    ax.plot(x_fine, g_inner_r, '-.', color='gray', lw=1, alpha=0.5)
    ax.plot(x_fine, g_center,  '-', color='blue', lw=1, alpha=0.4,
            label='center')
    ax.plot(x_fine, g_shoulder, '-', color='red', lw=1, alpha=0.4,
            label=f'shoulder ({mu4:.0f} cm)')

    ax.set_title(f'{temp} C', fontweight='bold')
    ax.set_xlim(0, 60)
    if idx == 0:
        ax.legend(fontsize=6, loc='upper left')
    if idx >= 4:
        ax.set_xlabel('Distance (cm)')
    ax.set_ylabel('T (C)')

for idx in range(len(nominal_temperatures), len(axes_flat)):
    axes_flat[idx].set_visible(False)

plt.suptitle('Global fit (ODR): individual profiles with components',
             fontweight='bold', fontsize=14)
plt.tight_layout()
plt.savefig('global_fit_odr_individual.png', dpi=150, bbox_inches='tight')
plt.show()

# %% [markdown]
# ## Comparison: LS vs ODR parameters

# %%
print("Comparison of parameters:")
print(f"{'param':>10s}  {'LS':>12s}  {'ODR':>12s}")
print("-" * 40)
param_names = ['baseline', 'a1', 'b1', 'mu1', 'sigma1',
               'a2', 'b2', 'mu2', 'sigma2',
               'a3', 'b3', 'sigma3',
               'a4', 'b4', 'mu4', 'sigma4']
for name, ls_val, odr_val in zip(param_names, best_p0, odr_result.beta):
    print(f"{name:>10s}  {ls_val:>12.4f}  {odr_val:>12.4f}")

# %% [markdown]
# ## Print the fitted function for use in other scripts

# %%
print("\ndef T_profile(x, nominal):")
print(f"    mid = {MIDPOINT}")
print(f"    return ({baseline:.4f}")
print(f"        + {a1:.6f} * (nominal - ({b1:.4f})) * np.exp(-0.5*((x - (mid - {mu1:.4f}))/{sigma1:.4f})**2)")
print(f"        + {a1:.6f} * (nominal - ({b1:.4f})) * np.exp(-0.5*((x - (mid + {mu1:.4f}))/{sigma1:.4f})**2)")
print(f"        + {a2:.6f} * (nominal - ({b2:.4f})) * np.exp(-0.5*((x - (mid - {mu2:.4f}))/{sigma2:.4f})**2)")
print(f"        + {a2:.6f} * (nominal - ({b2:.4f})) * np.exp(-0.5*((x - (mid + {mu2:.4f}))/{sigma2:.4f})**2)")
print(f"        + {a3:.6f} * (nominal - ({b3:.4f})) * np.exp(-0.5*((x - mid)/{sigma3:.4f})**2)")
print(f"        + {a4:.6f} * (nominal - ({b4:.4f})) * np.exp(-0.5*((x - {mu4:.4f})/{sigma4:.4f})**2))")

# %%
