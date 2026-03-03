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
# # Fitting Temperature Profiles with Sum of Gaussians
#
# Each measured temperature profile (from the literature) is modeled as
# the sum of 4 Gaussians centered near 15, 30, 30, and 45 cm,
# plus a constant baseline (room temperature).
#
# The optimizer finds the best heights, widths, and centers.

# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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

    # sort by distance
    order = np.argsort(distances)
    data[temp] = (distances[order], temperatures_at_distances[order])

# %% [markdown]
# ## Define the sum-of-Gaussians model
#
# 4 Gaussians (near 15, 30, 30, 45 cm) + a constant baseline.
# The two middle Gaussians at ~30 cm allow the model to capture
# the flat-topped peak seen in the profiles.

# %%
def sum_of_gaussians(x, baseline,
                     a1, mu1, sigma1,
                     a2, mu2, sigma2,
                     a3, mu3, sigma3,
                     a4, mu4, sigma4):
    """Sum of 4 Gaussians plus a constant baseline."""
    return (baseline
            + a1 * np.exp(-0.5 * ((x - mu1) / sigma1) ** 2)
            + a2 * np.exp(-0.5 * ((x - mu2) / sigma2) ** 2)
            + a3 * np.exp(-0.5 * ((x - mu3) / sigma3) ** 2)
            + a4 * np.exp(-0.5 * ((x - mu4) / sigma4) ** 2))

# %% [markdown]
# ## Fit each profile

# %%
fit_results = {}

for temp in nominal_temperatures:
    x_data, y_data = data[temp]
    peak_val = np.max(y_data)
    baseline_guess = 25.0

    # Initial guesses: baseline, then (amplitude, center, width) x 4
    # Gaussian 1: left shoulder ~15 cm
    # Gaussians 2,3: twin peaks at center ~27 and ~33 cm
    # Gaussian 4: right shoulder ~45 cm
    p0 = [baseline_guess,
           peak_val * 0.3, 15.0, 5.0,
           peak_val * 0.6, 27.0, 4.0,
           peak_val * 0.6, 33.0, 4.0,
           peak_val * 0.3, 45.0, 5.0]

    # Bounds: amplitudes > 0, centers within reactor, widths reasonable
    lower = [0,
             0, 8, 1.0,
             0, 20, 1.0,
             0, 26, 1.0,
             0, 38, 1.5]
    upper = [100,
             peak_val * 2, 22, 15,
             peak_val * 2, 33, 15,
             peak_val * 2, 40, 15,
             peak_val * 2, 52, 15]

    try:
        popt, pcov = curve_fit(sum_of_gaussians, x_data, y_data,
                               p0=p0, bounds=(lower, upper),
                               maxfev=50000)
        fit_results[temp] = popt
        residuals = y_data - sum_of_gaussians(x_data, *popt)
        rmse = np.sqrt(np.mean(residuals ** 2))
        print(f"T={temp} C  RMSE={rmse:.2f} C  "
              f"baseline={popt[0]:.1f}  "
              f"centers=[{popt[2]:.1f}, {popt[5]:.1f}, {popt[8]:.1f}, {popt[11]:.1f}]  "
              f"amplitudes=[{popt[1]:.1f}, {popt[4]:.1f}, {popt[7]:.1f}, {popt[10]:.1f}]  "
              f"widths=[{popt[3]:.1f}, {popt[6]:.1f}, {popt[9]:.1f}, {popt[12]:.1f}]")
    except RuntimeError as e:
        print(f"T={temp} C  FIT FAILED: {e}")
        fit_results[temp] = None

# %% [markdown]
# ## Plot: data points vs. fitted sum-of-Gaussians curves

# %%
colors = {800: '#4472C4',
          850: '#EE7D31',
          900: '#A5A5A5',
          950: '#FFC002',
          1000: '#5B9CD5',
          1050: '#70AE47',
          1100: '#264479'}

fig, ax = plt.subplots(figsize=(10, 6))

x_fine = np.linspace(0, 60, 300)

for temp in nominal_temperatures:
    x_data, y_data = data[temp]
    color = colors[temp]

    # data points
    ax.plot(x_data, y_data, 'o', color=color, markersize=4, alpha=0.8)

    # fitted curve
    if fit_results[temp] is not None:
        y_fit = sum_of_gaussians(x_fine, *fit_results[temp])
        ax.plot(x_fine, y_fit, '-', color=color, linewidth=1.5,
                label=f'{temp} C')

ax.set_xlabel('Distance into reactor (cm)', fontweight='bold')
ax.set_ylabel('Temperature (C)', fontweight='bold')
ax.set_title('Temperature profiles: data (points) vs. sum-of-Gaussians fit (lines)',
             fontweight='bold')
ax.legend(loc='upper left')
ax.set_xlim(0, 60)
ax.set_ylim(0, None)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontweight('bold')

plt.tight_layout()
plt.savefig('all_profiles_fit.png', dpi=150, bbox_inches='tight')
plt.show()

# %% [markdown]
# ## Individual profile plots with component Gaussians

# %%
fig, axes = plt.subplots(3, 3, figsize=(14, 12), sharex=True)
axes_flat = axes.flatten()

for idx, temp in enumerate(nominal_temperatures):
    ax = axes_flat[idx]
    x_data, y_data = data[temp]
    color = colors[temp]

    ax.plot(x_data, y_data, 'o', color='black', markersize=4, label='data')

    if fit_results[temp] is not None:
        popt = fit_results[temp]
        y_fit = sum_of_gaussians(x_fine, *popt)
        ax.plot(x_fine, y_fit, '-', color=color, linewidth=2, label='total fit')

        # plot individual Gaussians
        baseline = popt[0]
        gauss_labels = ['G1 (~15 cm)', 'G2 (~27 cm)', 'G3 (~33 cm)', 'G4 (~45 cm)']
        gauss_styles = ['--', '-.', ':', '--']
        for g in range(4):
            a = popt[1 + g * 3]
            mu = popt[2 + g * 3]
            sigma = popt[3 + g * 3]
            y_g = baseline + a * np.exp(-0.5 * ((x_fine - mu) / sigma) ** 2)
            ax.plot(x_fine, y_g, gauss_styles[g], color='gray',
                    linewidth=1, alpha=0.7, label=gauss_labels[g])

    ax.set_title(f'{temp} C', fontweight='bold')
    ax.set_xlim(0, 60)
    if idx == 0:
        ax.legend(fontsize=7, loc='upper left')
    if idx >= 4:
        ax.set_xlabel('Distance (cm)')
    ax.set_ylabel('T (C)')

# hide unused subplots
for idx in range(len(nominal_temperatures), len(axes_flat)):
    axes_flat[idx].set_visible(False)

plt.suptitle('Individual fits with component Gaussians', fontweight='bold', fontsize=14)
plt.tight_layout()
plt.savefig('individual_profile_fits.png', dpi=150, bbox_inches='tight')
plt.show()

# %% [markdown]
# ## Summary table of fit parameters

# %%
print(f"{'Temp':>6s}  {'baseline':>8s}  "
      f"{'a1':>7s} {'mu1':>5s} {'sig1':>5s}  "
      f"{'a2':>7s} {'mu2':>5s} {'sig2':>5s}  "
      f"{'a3':>7s} {'mu3':>5s} {'sig3':>5s}  "
      f"{'a4':>7s} {'mu4':>5s} {'sig4':>5s}")
print("-" * 110)

for temp in nominal_temperatures:
    if fit_results[temp] is not None:
        p = fit_results[temp]
        print(f"{temp:>6d}  {p[0]:>8.2f}  "
              f"{p[1]:>7.1f} {p[2]:>5.1f} {p[3]:>5.1f}  "
              f"{p[4]:>7.1f} {p[5]:>5.1f} {p[6]:>5.1f}  "
              f"{p[7]:>7.1f} {p[8]:>5.1f} {p[9]:>5.1f}  "
              f"{p[10]:>7.1f} {p[11]:>5.1f} {p[12]:>5.1f}")

# %%
