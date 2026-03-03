# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Global Fit of Temperature Profiles with Symmetric Sum of Gaussians
#
# All seven temperature profiles (800-1100 C) are fit simultaneously
# with a single set of 9 parameters. The profile is symmetric about
# the midpoint (30 cm), with Gaussian amplitudes that scale linearly
# with (nominal_temperature - offset).
#
# Model function:
#
# T(x, nominal) = baseline
#     + a1*(nominal - b1) * [G(x, midpoint-mu1, sigma1) + G(x, midpoint+mu1, sigma1)]
#     + a2*(nominal - b2) * [G(x, midpoint-mu2, sigma2) + G(x, midpoint+mu2, sigma2)]

# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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
# ## Define the symmetric sum-of-Gaussians model

# %%
MIDPOINT = 30.0

def sum_of_gaussians(x, nominal, baseline,
                     a1, b1, mu1, sigma1,
                     a2, b2, mu2, sigma2):
    """Symmetrical sum of 4 Gaussians plus a constant baseline.

    The 4 Gaussians are placed symmetrically about the midpoint (30 cm):
      - outer pair at midpoint +/- mu1, amplitude a1*(nominal - b1)
      - inner pair at midpoint +/- mu2, amplitude a2*(nominal - b2)
    """
    midpoint = MIDPOINT
    return (baseline
            + a1 * (nominal - b1) * np.exp(-0.5 * ((x - (midpoint - mu1)) / sigma1) ** 2)
            + a2 * (nominal - b2) * np.exp(-0.5 * ((x - (midpoint - mu2)) / sigma2) ** 2)
            + a2 * (nominal - b2) * np.exp(-0.5 * ((x - (midpoint + mu2)) / sigma2) ** 2)
            + a1 * (nominal - b1) * np.exp(-0.5 * ((x - (midpoint + mu1)) / sigma1) ** 2))

# %% [markdown]
# ## Build the global residual vector and optimize

# %%
# Stack all data for the global fit
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

def residuals(params):
    baseline, a1, b1, mu1, sigma1, a2, b2, mu2, sigma2 = params
    y_pred = sum_of_gaussians(x_all, nominal_all, baseline,
                               a1, b1, mu1, sigma1,
                               a2, b2, mu2, sigma2)
    return y_pred - y_all

# Multi-start optimization with physically motivated initial guesses.
#        baseline  a1         b1          mu1     sigma1   a2         b2          mu2     sigma2
bounds = [(0, 80), (0.001, 20), (-5000, 800), (8, 28), (1.5, 15), (0.001, 20), (-5000, 800), (0, 14), (1.5, 15)]

lower_b = [b[0] for b in bounds]
upper_b = [b[1] for b in bounds]

best_cost = np.inf
best_result = None

initial_guesses = [
    [25.0,  0.5,   100.0, 15.0, 4.0,  0.8,   100.0, 4.0, 4.5],
    [25.0,  0.3,  -100.0, 12.0, 3.0,  0.5,  -100.0, 3.0, 5.0],
    [30.0,  1.0,   300.0, 15.0, 5.0,  0.6,     0.0, 5.0, 4.0],
    [25.0,  0.2,  -200.0, 14.0, 4.5,  0.7,   200.0, 4.0, 5.0],
    [20.0,  0.8,   400.0, 12.0, 3.5,  0.4,  -300.0, 2.0, 6.0],
    [25.0,  0.15, -400.0, 13.0, 5.0,  0.5,     0.0, 3.5, 4.5],
    [30.0,  0.4,   200.0, 10.0, 5.5,  0.9,   300.0, 5.0, 3.5],
    [20.0,  2.0,   700.0, 10.0, 5.0,  0.3,  -800.0, 4.0, 4.0],
    [25.0,  0.1, -1000.0, 15.0, 5.0,  0.05,-1500.0, 4.0, 5.0],
    [25.0,  5.0,   750.0, 12.0, 4.0,  0.2, -1000.0, 3.0, 5.0],
    [25.0,  0.05,-3000.0, 12.0, 5.0,  0.03,-4000.0, 4.0, 4.0],
    [25.0,  3.0,   780.0, 10.0, 6.0,  0.1, -2000.0, 3.0, 4.0],
]

for i, p0 in enumerate(initial_guesses):
    try:
        res = least_squares(residuals, p0, bounds=(lower_b, upper_b),
                            method='trf', max_nfev=200000)
        if res.cost < best_cost:
            best_cost = res.cost
            best_result = res
            print(f"  Guess {i}: cost={res.cost:.1f} ** new best **")
        else:
            print(f"  Guess {i}: cost={res.cost:.1f}")
    except Exception as e:
        print(f"  Guess {i}: failed ({e})")

result = best_result
baseline, a1, b1, mu1, sigma1, a2, b2, mu2, sigma2 = result.x

print(f"\nBest optimization result (cost = {result.cost:.1f}):")
print(f"  baseline = {baseline:.2f} C")
print(f"  Outer Gaussians: a1={a1:.4f}, b1={b1:.1f}, mu1={mu1:.2f} cm, sigma1={sigma1:.2f} cm")
print(f"  Inner Gaussians: a2={a2:.4f}, b2={b2:.1f}, mu2={mu2:.2f} cm, sigma2={sigma2:.2f} cm")
print(f"\n  Gaussian centers (relative to midpoint={MIDPOINT} cm):")
print(f"    Outer pair at {MIDPOINT - mu1:.1f} and {MIDPOINT + mu1:.1f} cm")
print(f"    Inner pair at {MIDPOINT - mu2:.1f} and {MIDPOINT + mu2:.1f} cm")

# Per-profile RMSE
print()
for temp in nominal_temperatures:
    xd, yd = data[temp]
    y_pred = sum_of_gaussians(xd, temp, baseline, a1, b1, mu1, sigma1, a2, b2, mu2, sigma2)
    rmse = np.sqrt(np.mean((y_pred - yd) ** 2))
    print(f"  T={temp} C:  RMSE = {rmse:.1f} C")

# %% [markdown]
# ## Plot: all profiles, data (points) vs. global fit (lines)

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

    ax.plot(xd, yd, 'o', color=color, markersize=4, alpha=0.8)

    y_fit = sum_of_gaussians(x_fine, temp, baseline, a1, b1, mu1, sigma1, a2, b2, mu2, sigma2)
    ax.plot(x_fine, y_fit, '-', color=color, linewidth=1.5, label=f'{temp} C')

ax.set_xlabel('Distance into reactor (cm)', fontweight='bold')
ax.set_ylabel('Temperature (C)', fontweight='bold')
ax.set_title('Global fit: 9 shared parameters, symmetric about 30 cm', fontweight='bold')
ax.legend(loc='upper left')
ax.set_xlim(0, 60)
ax.set_ylim(0, None)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontweight('bold')

plt.tight_layout()
plt.savefig('global_fit_all_profiles.png', dpi=150, bbox_inches='tight')
plt.show()

# %% [markdown]
# ## Individual profile panels with component Gaussians

# %%
fig, axes = plt.subplots(3, 3, figsize=(14, 12), sharex=True)
axes_flat = axes.flatten()

for idx, temp in enumerate(nominal_temperatures):
    ax = axes_flat[idx]
    xd, yd = data[temp]
    color = colors[temp]

    ax.plot(xd, yd, 'o', color='black', markersize=4, label='data')

    y_fit = sum_of_gaussians(x_fine, temp, baseline, a1, b1, mu1, sigma1, a2, b2, mu2, sigma2)
    ax.plot(x_fine, y_fit, '-', color=color, linewidth=2, label='global fit')

    # Show the 4 component Gaussians (each shifted up by baseline for visual clarity)
    g_outer_left  = baseline + a1 * (temp - b1) * np.exp(-0.5 * ((x_fine - (MIDPOINT - mu1)) / sigma1) ** 2)
    g_inner_left  = baseline + a2 * (temp - b2) * np.exp(-0.5 * ((x_fine - (MIDPOINT - mu2)) / sigma2) ** 2)
    g_inner_right = baseline + a2 * (temp - b2) * np.exp(-0.5 * ((x_fine - (MIDPOINT + mu2)) / sigma2) ** 2)
    g_outer_right = baseline + a1 * (temp - b1) * np.exp(-0.5 * ((x_fine - (MIDPOINT + mu1)) / sigma1) ** 2)

    ax.plot(x_fine, g_outer_left,  '--', color='gray', linewidth=1, alpha=0.6, label=f'outer ({MIDPOINT-mu1:.0f} cm)')
    ax.plot(x_fine, g_inner_left,  '-.', color='gray', linewidth=1, alpha=0.6, label=f'inner ({MIDPOINT-mu2:.0f} cm)')
    ax.plot(x_fine, g_inner_right, '-.', color='gray', linewidth=1, alpha=0.6)
    ax.plot(x_fine, g_outer_right, '--', color='gray', linewidth=1, alpha=0.6)

    ax.set_title(f'{temp} C', fontweight='bold')
    ax.set_xlim(0, 60)
    if idx == 0:
        ax.legend(fontsize=7, loc='upper left')
    if idx >= 4:
        ax.set_xlabel('Distance (cm)')
    ax.set_ylabel('T (C)')

for idx in range(len(nominal_temperatures), len(axes_flat)):
    axes_flat[idx].set_visible(False)

plt.suptitle('Global fit: individual profiles with component Gaussians', fontweight='bold', fontsize=14)
plt.tight_layout()
plt.savefig('global_fit_individual.png', dpi=150, bbox_inches='tight')
plt.show()

# %% [markdown]
# ## Print the fitted function for use in other scripts

# %%
print("Fitted parameters:")
print(f"  MIDPOINT = {MIDPOINT}")
print(f"  baseline = {baseline:.4f}")
print(f"  a1 = {a1:.6f},  b1 = {b1:.4f},  mu1 = {mu1:.4f},  sigma1 = {sigma1:.4f}")
print(f"  a2 = {a2:.6f},  b2 = {b2:.4f},  mu2 = {mu2:.4f},  sigma2 = {sigma2:.4f}")
print()
print("def T_profile(x, nominal):")
print(f"    midpoint = {MIDPOINT}")
print(f"    return ({baseline:.4f}")
print(f"        + {a1:.6f} * (nominal - ({b1:.4f})) * np.exp(-0.5 * ((x - {MIDPOINT - mu1:.4f}) / {sigma1:.4f})**2)")
print(f"        + {a2:.6f} * (nominal - ({b2:.4f})) * np.exp(-0.5 * ((x - {MIDPOINT - mu2:.4f}) / {sigma2:.4f})**2)")
print(f"        + {a2:.6f} * (nominal - ({b2:.4f})) * np.exp(-0.5 * ((x - {MIDPOINT + mu2:.4f}) / {sigma2:.4f})**2)")
print(f"        + {a1:.6f} * (nominal - ({b1:.4f})) * np.exp(-0.5 * ((x - {MIDPOINT + mu1:.4f}) / {sigma1:.4f})**2))")

# %%
