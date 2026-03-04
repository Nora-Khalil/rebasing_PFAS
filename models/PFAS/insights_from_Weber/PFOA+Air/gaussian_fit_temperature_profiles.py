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
#     display_name: Python 3 (ipykernel)
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
import ipywidgets as widgets
from IPython.display import display
from scipy.odr import ODR, Model, RealData
from scipy.optimize import least_squares

# %%
# Load data
temperature_profiles = pd.read_csv('temperature_profiles_from_literature.csv')

nominal_temperatures = list(range(800, 1150, 50))
nominal_cooler_temperatures = list(range(400, 750, 50))

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
     a4, b4, mu4, sigma4
     ) = params

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

def report_parameters(params, cost=None):
    (baseline,
     a1, b1, mu1, sigma1,
     a2, b2, mu2, sigma2,
     a3, b3, sigma3,
     a4, b4, mu4, sigma4
     ) = params
    if cost:
       print(f"Cost: {cost:.4f}")
    print(f"\nFitted parameters:")
    print(f"  baseline = {baseline:.2f} C")
    print(f"  Outer pair:     a1={a1:.4f}, b1={b1:.1f}, mu1={mu1:.2f} cm, sigma1={sigma1:.2f} cm")
    print(f"                  centers at {MIDPOINT-mu1:.1f} and {MIDPOINT+mu1:.1f} cm")
    print(f"  Inner pair:     a2={a2:.4f}, b2={b2:.1f}, mu2={mu2:.2f} cm, sigma2={sigma2:.2f} cm")
    print(f"                  centers at {MIDPOINT-mu2:.1f} and {MIDPOINT+mu2:.1f} cm")
    print(f"  Center fill:    a3={a3:.4f}, b3={b3:.1f}, sigma3={sigma3:.2f} cm  (at {MIDPOINT} cm)")
    print(f"  Asym. shoulder: a4={a4:.4f}, b4={b4:.1f}, mu4={mu4:.2f} cm, sigma4={sigma4:.2f} cm")

# %%
colors = {800: '#4472C4',
          850: '#EE7D31',
          900: '#A5A5A5',
          950: '#FFC002',
          1000: '#5B9CD5',
          1050: '#70AE47',
          1100: '#264479'}

x_fine = np.linspace(0, 60, 300)

def plot_individual_profiles_with_components(params, errors=None):
    (baseline_, a1_, b1_, mu1_, sigma1_,
    a2_, b2_, mu2_, sigma2_,
    a3_, b3_, sigma3_,
    a4_, b4_, mu4_, sigma4_) = params
    fig, axes = plt.subplots(4, 4, figsize=(15, 12), sharex=True)
    axes_flat = axes.flatten()

    nominal_temperatures_to_plot = sorted(nominal_temperatures + nominal_cooler_temperatures, reverse=True)[:15]

    for idx, temp in enumerate(nominal_temperatures_to_plot):
        ax = axes_flat[idx]

        if temp in data:
            xd, yd = data[temp]
            color = colors[temp]

            if errors:
                sx_dist, sy_temp = errors

                ax.errorbar(xd, yd, xerr=sx_dist, yerr=sy_temp,
                            fmt='o', color='black', markersize=3, alpha=0.6,
                            elinewidth=0.5, capsize=0, label='data')
            else:
                ax.plot(xd, yd, 'o', color='black', markersize=3, alpha=0.6, label='data')

        y_fit = (baseline_
                 + a1_ * (temp - b1_) * np.exp(-0.5 * ((x_fine - (MIDPOINT - mu1_)) / sigma1_) ** 2)
                 + a1_ * (temp - b1_) * np.exp(-0.5 * ((x_fine - (MIDPOINT + mu1_)) / sigma1_) ** 2)
                 + a2_ * (temp - b2_) * np.exp(-0.5 * ((x_fine - (MIDPOINT - mu2_)) / sigma2_) ** 2)
                 + a2_ * (temp - b2_) * np.exp(-0.5 * ((x_fine - (MIDPOINT + mu2_)) / sigma2_) ** 2)
                 + a3_ * (temp - b3_) * np.exp(-0.5 * ((x_fine - MIDPOINT) / sigma3_) ** 2)
                 + a4_ * (temp - b4_) * np.exp(-0.5 * ((x_fine - mu4_) / sigma4_) ** 2))
        ax.plot(x_fine, y_fit, '-', color=color, linewidth=2, label='fit')

        b_off = baseline_
        g_outer_l = b_off + a1_ * (temp - b1_) * np.exp(-0.5 * ((x_fine - (MIDPOINT - mu1_)) / sigma1_) ** 2)
        g_outer_r = b_off + a1_ * (temp - b1_) * np.exp(-0.5 * ((x_fine - (MIDPOINT + mu1_)) / sigma1_) ** 2)
        g_inner_l = b_off + a2_ * (temp - b2_) * np.exp(-0.5 * ((x_fine - (MIDPOINT - mu2_)) / sigma2_) ** 2)
        g_inner_r = b_off + a2_ * (temp - b2_) * np.exp(-0.5 * ((x_fine - (MIDPOINT + mu2_)) / sigma2_) ** 2)
        g_center = b_off + a3_ * (temp - b3_) * np.exp(-0.5 * ((x_fine - MIDPOINT) / sigma3_) ** 2)
        g_shoulder = b_off + a4_ * (temp - b4_) * np.exp(-0.5 * ((x_fine - mu4_) / sigma4_) ** 2)

        ax.plot(x_fine, g_outer_l, '--', color='gray', lw=1, alpha=0.5,
                label=f'outer ({MIDPOINT-mu1_:.0f} cm)')
        ax.plot(x_fine, g_outer_r, '--', color='gray', lw=1, alpha=0.5)
        ax.plot(x_fine, g_inner_l, '-.', color='gray', lw=1, alpha=0.5,
                label=f'inner ({MIDPOINT-mu2_:.0f} cm)')
        ax.plot(x_fine, g_inner_r, '-.', color='gray', lw=1, alpha=0.5)
        ax.plot(x_fine, g_center, '-', color='blue', lw=1, alpha=0.4, label='center')
        ax.plot(x_fine, g_shoulder, '-', color='red', lw=1, alpha=0.4,
                label=f'shoulder ({mu4_:.0f} cm)')

        ax.set_title(f'{temp} C', fontweight='bold')
        ax.set_xlim(0, 60)
        ax.set_ylim(0, 1200)
        if idx == 0:
            ax.legend(fontsize=6, loc='upper left')
        if idx >= 4:
            ax.set_xlabel('Distance (cm)')
        ax.set_ylabel('T (C)')

    for idx in range(len(nominal_temperatures_to_plot), len(axes_flat)):
        axes_flat[idx].set_visible(False)

    plt.suptitle('Individual profiles with components',
                 fontweight='bold', fontsize=14)
    plt.tight_layout()
    plt.show()
plot_individual_profiles_with_components(best_p0)

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


# %%

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

lower = np.array([0, 0.001, 0,  5, 3.,
                     0.001, 0,  0, 3.,
                     0.001, 0,     3.,
                     0.001, 0,  0, 3.])
upper = np.array([80,   20,   1100, 30, 25, # outer pair
                        20,   1100, 15, 25, # inner pair
                        20,   1100,     25, # center fill
                        20,   1100, 30, 25]) # shoulder

best_cost = np.inf
best_p0 = None

initial_guesses = [
#baseline, a1, b1, mu1, sigma1, a2, b2, mu2, sigma2, a3, b3, sigma3, a4, b4, mu4, sigma4
]
# sobol sequence between lower and upper bounds
from scipy.stats import qmc
sampler = qmc.Sobol(d=len(lower), scramble=True)
n_samples = 20
sobol_samples = sampler.random(n_samples)
initial_guesses = qmc.scale(sobol_samples, lower, upper)
initial_guesses = initial_guesses.round(decimals=1)
# display as decimal not float
initial_guesses[:, 0] = 25.0
pd.DataFrame(initial_guesses, columns=['baseline', 'a1', 'b1', 'mu1', 'sigma1', 'a2', 'b2', 'mu2', 'sigma2', 'a3', 'b3', 'sigma3', 'a4', 'b4', 'mu4', 'sigma4'])



# %%
print("Finding good starting parameters via least squares...")
best_cost = np.inf
for i, p0 in enumerate(initial_guesses):
    try:
        res = least_squares(residuals_ls, p0, bounds=(lower, upper),
                            method='trf', max_nfev=50000)
        if res.cost < best_cost:
            best_cost = res.cost
            best_p0 = res.x.copy()
            print(f"  Guess {i}: cost={res.cost:.1f} ** new best **")
            plot_individual_profiles_with_components(res.x)
        else:
            print(f"  Guess {i}: cost={res.cost:.1f}")
    except Exception as e:
        print(f"  Guess {i}: failed ({e})")

print(f"\nBest LS cost = {best_cost:.1f}")
print(f"Starting ODR from these parameters.\n")

# %%
report_parameters(best_p0, cost=best_cost)

# %%
plot_individual_profiles_with_components(best_p0)

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

params = odr_result.beta

report_parameters(params)

# Per-profile RMSE (in T only, for comparison with earlier fits)
print()
for temp in nominal_temperatures:
    xd, yd = data[temp]
    xn = np.vstack([xd, np.full_like(xd, temp)])
    y_pred = sum_of_gaussians_odr(odr_result.beta, xn)
    rmse = np.sqrt(np.mean((y_pred - yd) ** 2))
    print(f"  T={temp} C:  RMSE = {rmse:.1f} C")

# %%
report_parameters(odr_result.beta, cost=odr_result.res_var)

# uncertainties
report_parameters(odr_result.sd_beta)

# %%
plot_individual_profiles_with_components(odr_result.beta, errors=(sx_dist, sy_temp))

# %% [markdown]
# ## Plot: all profiles overlay

# %%


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

for temp in nominal_cooler_temperatures:
    # fit line only
    xn_fine = np.vstack([x_fine, np.full_like(x_fine, temp)])
    y_fit = sum_of_gaussians_odr(odr_result.beta, xn_fine)
    ax.plot(x_fine, y_fit, '--', linewidth=1.0, label=f'{temp} C')

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

# %%


slider_style = {'description_width': '80px'}
layout = widgets.Layout(width='380px')

w_baseline = widgets.FloatSlider(value=float(baseline), min=0.0, max=120.0, step=0.1, description='baseline', style=slider_style, layout=layout)
w_a1 = widgets.FloatSlider(value=float(a1), min=0.001, max=20.0, step=0.001, description='a1', style=slider_style, layout=layout)
w_b1 = widgets.FloatSlider(value=float(b1), min=-5000.0, max=800.0, step=1.0, description='b1', style=slider_style, layout=layout)
w_mu1 = widgets.FloatSlider(value=float(mu1), min=5.0, max=28.0, step=0.1, description='mu1', style=slider_style, layout=layout)
w_sigma1 = widgets.FloatSlider(value=float(sigma1), min=1.5, max=15.0, step=0.1, description='sigma1', style=slider_style, layout=layout)

w_a2 = widgets.FloatSlider(value=float(a2), min=0.001, max=20.0, step=0.001, description='a2', style=slider_style, layout=layout)
w_b2 = widgets.FloatSlider(value=float(b2), min=-5000.0, max=800.0, step=1.0, description='b2', style=slider_style, layout=layout)
w_mu2 = widgets.FloatSlider(value=float(mu2), min=0.0, max=14.0, step=0.1, description='mu2', style=slider_style, layout=layout)
w_sigma2 = widgets.FloatSlider(value=float(sigma2), min=1.5, max=15.0, step=0.1, description='sigma2', style=slider_style, layout=layout)

w_a3 = widgets.FloatSlider(value=float(a3), min=0.001, max=20.0, step=0.001, description='a3', style=slider_style, layout=layout)
w_b3 = widgets.FloatSlider(value=float(b3), min=-5000.0, max=800.0, step=1.0, description='b3', style=slider_style, layout=layout)
w_sigma3 = widgets.FloatSlider(value=float(sigma3), min=1.5, max=15.0, step=0.1, description='sigma3', style=slider_style, layout=layout)

w_a4 = widgets.FloatSlider(value=float(a4), min=0.001, max=20.0, step=0.001, description='a4', style=slider_style, layout=layout)
w_b4 = widgets.FloatSlider(value=float(b4), min=-5000.0, max=800.0, step=1.0, description='b4', style=slider_style, layout=layout)
w_mu4 = widgets.FloatSlider(value=float(mu4), min=8.0, max=52.0, step=0.1, description='mu4', style=slider_style, layout=layout)
w_sigma4 = widgets.FloatSlider(value=float(sigma4), min=1.5, max=15.0, step=0.1, description='sigma4', style=slider_style, layout=layout)

controls = widgets.HBox([
    widgets.VBox([w_baseline, w_a1, w_b1, w_mu1, w_sigma1, w_a2, w_b2, w_mu2]),
    widgets.VBox([w_sigma2, w_a3, w_b3, w_sigma3, w_a4, w_b4, w_mu4, w_sigma4])
])

out = widgets.interactive_output(
    plot_individual_profiles_with_components,
    {
        'baseline_': w_baseline, 'a1_': w_a1, 'b1_': w_b1, 'mu1_': w_mu1, 'sigma1_': w_sigma1,
        'a2_': w_a2, 'b2_': w_b2, 'mu2_': w_mu2, 'sigma2_': w_sigma2,
        'a3_': w_a3, 'b3_': w_b3, 'sigma3_': w_sigma3,
        'a4_': w_a4, 'b4_': w_b4, 'mu4_': w_mu4, 'sigma4_': w_sigma4
    }
)

display(widgets.HTML('<h4>Adjust fit parameters</h4>'))
display(controls, out)

# %% [markdown]
# ## Comparison: LS vs ODR parameters

# %%

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
