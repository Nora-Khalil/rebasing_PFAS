# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: ct32
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Plug flow reactor simulation for PFOA thermal degradation
#
# Uses a Lagrangian particle approach with an ExtensibleIdealGasConstPressureMoleReactor.
# The temperature profile T(x) is imposed via the extensible reactor API, overriding the
# energy equation RHS. Position x is added as an extra state variable with dx/dt = u.
#
# The temperature profile is a sum of Gaussians fitted to the experimental data from
# Weber et al., parameterized by the nominal furnace temperature.
#
#

# %%
import os
import re
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# %% [markdown]
#
# Temperature profile: T(x_cm, nominal_C) -> T in °C
# Fitted sum of 5 Gaussians (ODR fit from gaussian_fit_temperature_profiles.py)
# x in cm, nominal in °C, returns °C
#
#

# %%
# Parameters from the ODR fit
_BASELINE = 21.9537   # °C
_A1 = 1.942512;  _B1 = 925.0228;  _MU1 = 11.9773;  _SIGMA1 = 3.9644
_A2 = 0.630055;  _B2 = 0.0000;    _MU2 = 4.3965;   _SIGMA2 = 4.0971
_A3 = 0.343392;  _B3 = 447.7876;  _SIGMA3 = 10.7594
_MIDPOINT = 30.0  # cm


# %%
def T_profile_celsius(x_cm, nominal_C):
    """Temperature profile in °C, given position in cm and nominal temp in °C."""
    x = x_cm
    nom = nominal_C
    mid = _MIDPOINT
    return (_BASELINE
        + _A1 * max(nom - _B1, 0) * np.exp(-0.5 * ((x - (mid - _MU1)) / _SIGMA1)**2)
        + _A1 * max(nom - _B1, 0) * np.exp(-0.5 * ((x - (mid + _MU1)) / _SIGMA1)**2)
        + _A2 * max(nom - _B2, 0) * np.exp(-0.5 * ((x - (mid - _MU2)) / _SIGMA2)**2)
        + _A2 * max(nom - _B2, 0) * np.exp(-0.5 * ((x - (mid + _MU2)) / _SIGMA2)**2)
        + _A3 * max(nom - _B3, 0) * np.exp(-0.5 * ((x - mid) / _SIGMA3)**2)
    )


# %%
def T_profile_K(x_m, nominal_C):
    """Temperature in K, given position in m and nominal temp in °C."""
    x_cm = x_m * 100.0  # m -> cm
    return T_profile_celsius(x_cm, nominal_C) + 273.15

# %%
def T_average_K(nominal_C):
    """Temperature in K averaged over the reactor length, given nominal temp in °C.
    
    But weight it differently, so that it gets the RTD better.
    do 1/(mean of 1/T) instead of mean of T
    """

    return 1/(np.mean([1/T_profile_K(x, nominal_C) for x in np.linspace(0, _MIDPOINT*2, 100)]))



# %%
def dTdx_K_per_m(x_m, nominal_C):
    """Analytical dT/dx in K/m, given position in m and nominal temp in °C.

    Each Gaussian term a * coeff * exp(-0.5*((x_cm - c)/s)^2) contributes
    dT/dx_cm = a * coeff * [-(x_cm - c)/s^2] * exp(-0.5*((x_cm - c)/s)^2)

    Then dT/dx_m = dT/dx_cm * 100 (chain rule: dx_cm/dx_m = 100)
    and the °C -> K conversion doesn't affect the derivative.
    """
    x_cm = x_m * 100.0
    nom = nominal_C
    mid = _MIDPOINT

    dTdx_cm = 0.0

    # Helper: derivative of a * coeff * exp(-0.5*((x - c)/s)^2) w.r.t. x
    def gauss_deriv(a, coeff, c, s):
        z = (x_cm - c) / s
        return a * coeff * (-z / s) * np.exp(-0.5 * z**2)

    coeff1 = max(nom - _B1, 0)
    coeff2 = max(nom - _B2, 0)
    coeff3 = max(nom - _B3, 0)

    # Outer pair
    dTdx_cm += gauss_deriv(_A1, coeff1, mid - _MU1, _SIGMA1)
    dTdx_cm += gauss_deriv(_A1, coeff1, mid + _MU1, _SIGMA1)
    # Inner pair
    dTdx_cm += gauss_deriv(_A2, coeff2, mid - _MU2, _SIGMA2)
    dTdx_cm += gauss_deriv(_A2, coeff2, mid + _MU2, _SIGMA2)
    # Center
    dTdx_cm += gauss_deriv(_A3, coeff3, mid, _SIGMA3)

    # Convert cm -> m: dT/dx_m = dT/dx_cm * (dx_cm / dx_m) = dT/dx_cm * 100
    return dTdx_cm * 100.0


# %% [markdown]
#
# ## Extensible Reactor: 
# imposed T(x) profile with x as a state variable
#

# %%
class PFR_Reactor(ct.ExtensibleIdealGasConstPressureMoleReactor):
    """Lagrangian particle reactor with imposed temperature profile T(x).

    Extra state variable: x (position in m).
    The energy equation RHS is overwritten to impose dT/dt = dT/dx * u.
    An additional equation dx/dt = u is added.
    """

    def __init__(self, *args, nominal_temp_C, mass_flow_rate, cross_section_area,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.nominal_temp_C = nominal_temp_C
        self.mass_flow_rate = mass_flow_rate      # kg/s (constant)
        self.cross_section_area = cross_section_area  # m^2 (constant)
        self.position = 0.0  # current position in m

    def after_initialize(self, t0):
        """Called after initialize. Register extra state variable."""
        if t0 != 0:
            raise ValueError(f"Expected initialization time t0=0, got t0={t0}")
        # Add one extra equation for position x
        self.n_vars += 1
        self.i_position = self.n_vars - 1  # index of the position variable
        self.position = 0.0

    def after_get_state(self, y):
        """Append x to the state vector."""
        y[self.n_vars - 1] = self.position

    def after_update_state(self, y):
        """Read x back from the state vector."""
        self.position = y[self.i_position]

    def before_component_index(self, name):
        # Other components are handled by the method from the base Reactor class
        if name == 'position':
            return self.i_position

    def before_component_name(self, i):
        # Other components are handled by the method from the base Reactor class
        if i == self.i_position:
            return 'position'


    def after_eval(self, t, LHS, RHS):
        """Modify the RHS to impose the temperature profile and add dx/dt = u."""
        # Current position from the state
        x_m = self.position

        # Current velocity: u = mdot / (rho * A)
        rho = self.density
        u = self.mass_flow_rate / (rho * self.cross_section_area)

        # Impose dT/dt = dT/dx * u
        dTdx = dTdx_K_per_m(x_m, self.nominal_temp_C)
        dTdt = dTdx * u

        # Temperature is state variable index 0 for IdealGasConstPressureMoleReactor
        # The governing eq is: LHS[0] * dT/dt = RHS[0]
        # So to impose our dT/dt: RHS[0] = LHS[0] * dT/dt
        RHS[0] = LHS[0] * dTdt

        # Extra equation for position: dx/dt = u
        # Extra variables start after the base state variables.
        # LHS for extra variables defaults to 1.0, so RHS = dx/dt = u
        RHS[self.n_vars - 1] = u


# %% [markdown]
#
# ## Simulation parameters
# (Weber et al. experiments)
#

# %%
# Jupyter-safe: __file__ is not defined in notebooks
here = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
mechanism_file = os.path.abspath(
    os.path.join(here, '..', 'chemkin_for_Weber_2025_12_12', 'PFAS_O2_H2O_N2O_publish.yaml')
)
print(mechanism_file)
assert os.path.isfile(mechanism_file), f"Mechanism file not found: {mechanism_file}"

# %%
# Reactor geometry
int_diam = 0.007       # m (7 mm ID)
radius = int_diam / 2
cross_area = np.pi * radius * radius  # m^2

# %%
# Flow conditions
Vdot = 150e-6 / 60.0  # volumetric flow rate: 150 mL/min -> m^3/s (at inlet T)

# %%
# Reactor length
length = 0.60  # m (60 cm)

# %%
# Initial composition (mole fractions)
initial_composition = {
    "PFOA": 4.02e-4,   # 402 ppm
    "H2O": 7.50e-4,    # 750 ppm
    "O2": 0.20975808,
    "N2": 0.78908992,
}

# %%
# Nominal furnace temperatures to simulate (°C)
nominal_temperatures = [400, 500, 600, 700, 800, 900, 1000, 1100]

# %%
# Number of time steps for output
n_steps = 200


# %%

def plot_rates_of_progress(gas):
    """Plot forward (blue) and reverse (red) rates of progress for all reactions."""
    # Plot forward rates
    plt.plot(gas.forward_rates_of_progress, 'b.', label='Forward')

    # Plot reverse rates
    plt.plot(gas.reverse_rates_of_progress, 'r.', label='Reverse')

    # label the highest of each
    top_forward = np.argsort(gas.forward_rates_of_progress)[-3:]
    top_reverse = np.argsort(gas.reverse_rates_of_progress)[-3:]
    for idx in top_forward:
        plt.text(idx, gas.forward_rates_of_progress[idx], f'R{idx}', color='blue', fontsize=8)
        print(f"R{idx}: {gas.reaction(idx).equation}  Forward ROP: {gas.forward_rates_of_progress[idx]:.3e}")
    for idx in top_reverse:
        plt.text(idx, gas.reverse_rates_of_progress[idx], f'R{idx}', color='red', fontsize=8)
        print(f"R{idx}: {gas.reaction(idx).equation}  Reverse ROP: {gas.reverse_rates_of_progress[idx]:.3e}")

    plt.xlabel('Reaction')
    plt.ylabel('Rate of Progress')
    plt.title('Rates of Progress for All Reactions')
    plt.legend()
    plt.show()

# %%
def find_culprit_reactions(gas, species_index=None, species_name=None):
    """Find the top reactions contributing to the net production/consumption of a species."""
    if species_index is not None:
        k = species_index
    elif species_name is not None:
        k = gas.species_index(species_name)
    else:
        raise ValueError("Either species_index or species_name must be provided")

    # Stoich coefficients for species k across all reactions
    nu = gas.product_stoich_coeffs[k, :] - gas.reactant_stoich_coeffs[k, :]
    
    # Absolute total rate
    rop = abs(gas.forward_rates_of_progress) + abs(gas.reverse_rates_of_progress)
    
    # Contribution of each reaction to species k
    contrib = nu * rop
    
    # Sort by magnitude
    top = np.argsort(np.abs(contrib))[::-1][:8]
    
    print(f"\nTop reactions affecting {gas.species_name(k)}:")
    for idx in top:
        if contrib[idx] != 0:
            print(f"  R{idx}: {gas.reaction(idx).equation}  "
                  f"(rate contribution: {contrib[idx]:.4e})")


# %% [markdown]
#
# # Run simulations

# %%
results = {}
gas = ct.Solution(mechanism_file)

# %%
# Mass flow rate (constant throughout the reactor)
# Vdot = 150 mL/min is at room temperature (before entering furnace).
# So compute mass flow rate at ~298 K, 1 atm:
gas.TPX = 298.15, ct.one_atm, initial_composition
rho_inlet = gas.density
mass_flow_rate = rho_inlet * Vdot


for nominal_T_C in nominal_temperatures:
    print(f"\n{'='*60}")
    print(f"Simulating nominal temperature: {nominal_T_C} °C")
    print(f"{'='*60}")

    # Set up the gas at inlet conditions
    T_inlet_K = T_profile_K(0.0, nominal_T_C)
    gas.TPX = T_inlet_K, ct.one_atm, initial_composition

    print(f"  Inlet T: {T_inlet_K:.1f} K ({T_inlet_K - 273.15:.1f} °C)")
    print(f"  Mass flow rate: {mass_flow_rate:.6e} kg/s")
    print(f"  Inlet velocity (at 298.15K): {Vdot / cross_area:.4f} m/s")
    print(f"  Inlet velocity (at T_inlet): {mass_flow_rate / (gas.density * cross_area):.4f} m/s")

    # Create the extensible reactor
    reactor = PFR_Reactor(
        gas, clone=True,
        nominal_temp_C=nominal_T_C,
        mass_flow_rate=mass_flow_rate,
        cross_section_area=cross_area,
    )

    # Create reactor network
    sim = ct.ReactorNet([reactor])

    # Estimate total residence time for setting up time steps
    # Use a rough average velocity
    # We do a distance-weighted average but should have a time-weighted
    T_avg_K = T_average_K(nominal_T_C)
    rho_avg = gas.density * T_inlet_K / T_avg_K  # ideal gas approximation
    u_avg = mass_flow_rate / (rho_avg * cross_area)
    t_total_est = length / u_avg
    dt = t_total_est / n_steps

    print(f"  Estimated avg velocity: {u_avg:.4f} m/s")
    print(f"  Estimated residence time: {t_total_est:.6f} s")
    print(f"  Time step: {dt:.6e} s")

    # Storage for results
    t_array = np.zeros(n_steps)
    x_array = np.zeros(n_steps)
    T_array = np.zeros(n_steps)
    states = ct.SolutionArray(reactor.phase)

    # Integrate
    for n in range(n_steps):
        t_i = (n + 1) * dt
        try:
            sim.advance(t_i)
        except ct.CanteraError as e:
            print(f"  Error occurred at step {n+1}: {e}")
            m = re.search(r"Components with largest weighted error estimates:\n(\d+).*\n(\d+).*\n(\d+)", str(e))
            if m:
                print(f"  Problematic components: {m.group(1)}, {m.group(2)}, {m.group(3)}")
                print(sim.component_name(int(m.group(1))))
                print(sim.component_name(int(m.group(2))))
                print(sim.component_name(int(m.group(3))))
                t_array = t_array[:n] # truncate
                x_array = x_array[:n]
                T_array = T_array[:n]
                last_good_state = states[-1]
                find_culprit_reactions(last_good_state, species_index=int(m.group(1))-1)
                find_culprit_reactions(last_good_state, species_index=int(m.group(2))-1)
                find_culprit_reactions(last_good_state, species_index=int(m.group(3))-1)
                plot_rates_of_progress(last_good_state)
            break

        t_array[n] = sim.time
        x_array[n] = reactor.position
        T_array[n] = reactor.phase.T
        states.append(reactor.phase.state)

        # Stop if we've passed the end of the reactor
        if reactor.position >= length:
            t_array = t_array[:n+1]
            x_array = x_array[:n+1]
            T_array = T_array[:n+1]
            print(f"  Reached end of reactor at step {n+1}, "
                  f"t = {sim.time:.6f} s, x = {reactor.position:.4f} m")
            break
    else:
        print(f"  Final position: x = {reactor.position:.4f} m (target: {length:.4f} m)")
        if reactor.position < length * 0.95:
            print(
                "  WARNING: Did not reach end of reactor. "
                "Consider increasing n_steps or t_total_est."
            )

    results[nominal_T_C] = {
        't': t_array,
        'x': x_array,
        'T': T_array,
        'states': states,
    }


# %% [markdown]
#
# # Plot results
#

# %%
colors = {
    400: '#4472C4', 500: '#EE7D31', 600: '#A5A5A5', 700: '#FFC002',
    800: '#5B9CD5', 900: '#70AE47', 1000: '#264479', 1100: '#9E480E',
    }

# %%
# --- Temperature profiles: imposed vs simulated ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
x_fine = np.linspace(0, length, 500)
for nominal_T_C in nominal_temperatures:
    color = colors[nominal_T_C]
    # Analytical profile
    T_analytical = np.array([T_profile_K(x, nominal_T_C) for x in x_fine]) - 273.15
    ax1.plot(x_fine * 100, T_analytical, '-', color=color, linewidth=1, alpha=0.5)
    # Simulated
    r = results[nominal_T_C]
    ax1.plot(r['x'] * 100, r['T'] - 273.15, '--', color=color, linewidth=1.5,
             label=f'{nominal_T_C} °C')
ax1.set_xlabel('Distance (cm)')
ax1.set_ylabel('Temperature (°C)')
ax1.set_title('Temperature profiles: target (solid) vs simulated (dashed)')
ax1.legend()
ax1.set_xlim(0, 60)
# --- PFOA destruction ---
for nominal_T_C in nominal_temperatures:
    color = colors[nominal_T_C]
    r = results[nominal_T_C]
    states = r['states']
    if 'PFOA' in states.species_names:
        idx = states.species_index('PFOA')
        X_PFOA = states.X[:, idx]
        X_PFOA_0 = X_PFOA[0] if X_PFOA[0] > 0 else 4.02e-4
        ax2.plot(r['x'] * 100, X_PFOA / X_PFOA_0, '-', color=color,
                 label=f'{nominal_T_C} °C')

ax2.set_xlabel('Distance (cm)')
ax2.set_ylabel('[PFOA] / [PFOA]₀')
ax2.set_title('PFOA along reactor')
ax2.legend()
ax2.set_xlim(0, 60)
ax2.set_ylim(0, 1.1)

plt.tight_layout()
plt.show()
plt.savefig('pfr_lagrangian_results.png', dpi=300)
# %%
# HF concentration along the reactor
plt.figure(figsize=(5, 4))
ax = plt.gca()
for nominal_T_C in nominal_temperatures:
    color = colors[nominal_T_C]
    r = results[nominal_T_C]
    states = r['states']
    if 'HF' in states.species_names:
        idx = states.species_index('HF')
        X_HF = states.X[:, idx]
        ax.plot(r['x'] * 100, X_HF, '-', color=color,
                 label=f'{nominal_T_C} °C')
ax.set_xlabel('Distance (cm)')
ax.set_ylabel('Mole fraction of HF')
ax.set_title('HF formation along reactor')
ax.legend()
ax.set_xlim(0, 60)
plt.tight_layout()
plt.show()

# %%
