# Exo-FMS_column_nongrey

Major Update History:
 - May 2021 - initial models
 - Dec 2021 - major overhaul
 - May 2022 - SW zenith angle geometric correction added (not ready for all modes yet)
 - Jun 2022 - Bezier short char added
 - Aug 2023 - VIM added - removed some methods + refactored arrays - much more stable

Elspeth K.H. Lee - Aug 2023

This is one part of a series of codes that build upon different two-stream approaches and schemes, primarily useful for the GCM modelling community.
This is the non-grey, picket fence version, with three bands in the visible, representing incident radiation from a host star, and two bands in the IR, representing the internal radiation propagating inside the planetary atmosphere.
This is for gas giant planets only with a defined internal temperature, a solid/liquid surface RT is not included here (but could be added easily if needed for special cases).

Here we colloquially use "two-stream" as the general methodology, making specific wether a particular method uses two-stream or four-stream (or n-stream) for it's calculations.

Some useful references useful for non-grey, picket fence modelling are: \
Chandrasekhar (1935) \
Parmentier & Guillot (2014) \
Parmentier et al.  (2015) \
Lee et al. (2021)

The non-grey, picket fence scheme builds upon the semi-grey scheme, by having 3 bands in the visible, and 2 bands in the IR.
The 3 V bands attempt to better represent the deposition of stellar flux into the atmosphere, and therefore produce a more realistic stellar heating profile.
The 2 IR bands (picket fence scheme), are used to both represent the line opacity in one band, and the continuum opacity in the other.
The fraction of the total flux in each band is controlled by a parameter (Beta_IR).
This allows the atmosphere to cool from two photospheric regions, rather than one in the semi-grey scheme, leading to more realistic cooling in the upper atmosphere.
However, at very low pressures, isothermal T-p profiles are still produced.

Using the fitting functions and relations from Parmentier et al. (2014, 2015) a quite realistic T-p profile can be produced.
This is because the fitting functions were designed to reproduce the results of a correlated-k scheme as best as possible.
Unfortunately the numerical results in 1D can in certain circumstances produce small spiky profiles, primarily due to the sensitivity of the opacity fitting function to the temperature, or, when the fitting function is outside it's valid range of temperatures (~100-4000 K) or pressures (~1e-4-1e3 pa).
Inside a GCM these small spikes are typically smoothed out by the dynamical processes.

To compile enter 'make' in the main directory. To remove compiled code enter 'make clean'. \
Some compiler options for gfortran, nvfortran and ifort are provided in the makefile.

This code performs various two-stream approaches from the literature in a non-grey, picket fence context:
1. Isothermal layer approximation
2. Toon method (with scattering and without scattering versions)
3. Short Characteristics method (linear or Bezier interpolant versions)
4. Heng et al. method (Improved Heng TS method in development)
5. Two-stream DISORT version (with T_int modifications by Xianyu Tan)
6. Variational Iteration Method (VIM) with analytical LW scattering

You can also see the header comments in the source code for some additional information.

For the shortwave fluxes, for methods that do not contain a shortwave scattering mode we include the 'adding method' (Mendonca et al. 2015 + references).
We detect if any albedo is present in the column, and perform the adding method to calculate the scattered flux, otherwise if there is zero albedo only the direct beam is used.

This emulates a single column inside the Exo-FMS GCM and is useful for testing and developing new techniques
as they would perform inside a GCM setting. This is also useful to see differences in each method and their various approximations.

We also include dry convective adjustment schemes, currently only 'Ray_adj', based on Raymond Pierrehumbert's python code.

# Namelist options

In the file 'FMS_RC.nml' you can select different options that control the simulation

ts_scheme: \
'Isothermal' - Isothermal ts method \
'Isothermal_2' - Isothermal ts method - high optical depth version \
'Toon' - Toon et al. method - four stream \
'Toon_scatter' - Toon et al. ts method with scattering - four stream \
'Shortchar_linear' - Short characteristics method with linear interpolants - four stream \
'Shortchar_Bezier' - Short characteristics method with Bezier interpolants - four stream  \
'Heng' - Heng et al. method \
'Disort_scatter' - two-stream DISORT version with scattering
'VIM' - Variational Iteration Method (VIM) - four stream

opac_scheme: \
'Freedman' - Uses the Rosseland mean fitting function from Freedman et al. (2014) \
'Valencia' - Uses the Rosseland mean fitting function from Valencia et al. (2013) \

adj_scheme: \
'Ray_dry' - Ray Pierrehumbert's dry convective adjustment scheme

The option 'None' for each of these scheme will it off (e.g. To run without conv adjustment set adj_scheme = 'None')

nlay, a_sh , b_sh - the number of layers, and filenames that contain the a and b constants for the hybrid sigma grid

pref - reference surface pressure (pa)

t_step - time step in seconds \
nstep - number of integer timesteps \
Rd_air - specific gas constant of the air (J kg-1 K-1) \
cp_air - specific heat capacity (constant pressure) of the air (J kg-1 K-1) \
grav - gravitational acceleration constant (m s-2) \
mu_z - cosine angle of the solar zenith angle \
Tirr - Irradiation temperature \
Tint - Internal temperature

k_V - visible band opacity (m2 kg-1) \
k_IR - IR band opacity (m2 kg-1) \
AB - Bond albedo \
fl - The Heng et al. (2011) parameter used for pressure dependent IR optical depths \
met - metallicity in dex solar (M/H)

Bezier - use Bezier interpolation for Temperature layer to level interpolation (.True.)

sw_ac(3) - shortwave single scattering albedo (constant at all layers) \
sw_gc(3) - shortwave asymmetry factor (constant at all layers) \
lw_ac(2) - longwave single scattering albedo (constant at all layers) \
lw_gc(2) - longwave asymmetry factor (constant at all layers)

zcorr - include zenith angle correction (.True.) \
zcorr_meth - zenith angle correction method (1,2)  \
radius - radius of the planet at surface (m)

iIC - Initial condition selection integer (4 = Parmentier et al. (2015) profile using Tirr and Tint) \
corr - Flag to perform the adiabatic gradient correction in the initial conditions \
table_num - 1 = Parmentier et al. (2015) table with TiO & VO, 2 = without TiO & VO table

# Plotting results

A python plotting routine, 'plot_TP.py' is provided to plot the results of the simulation.

# Gaussian Ordinance values

In methods that use Gaussian quadrature to perform the mu integration in intensity (to calculate the flux), various values and weights can be changed for testing at the top of the two-stream modules.
You can use these to switch between two-stream, four-stream etc etc
You will need to clean and recompile the code if these are changed.

# Personal recommendations

For non-scattering problems, we generally recommend that the short characteristics method be used with linear, as it is fast, efficient, very stable and also very accurate. This is currently what is used inside Exo-FMS for the Hot Jupiter simulations, and is even fast enough for high spatial resolution cases.
For shortwave scattering problems we recommend the adding method as included (or using the four-stream Toon or two-stream DISORT methods), the adding method is generally fast and accurate (enough), especially for grey opacity problems.
For longwave scattering problems we recommend the VIM code, as a fast but accurate method due to it's analytical properties.
If that fails try the Toon_scattering or DISORT code, they are generally reliable but slower compared to other scattering methods.

# Future developments

Additional opacity schemes from the literature can be added. \
Interpolating directly from the Freedman et al. (2014) tables rather than using the fitting function. \
Adding analytical four stream shortwave methods \
Adding more analytical longwave scattering methods
