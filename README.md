# Exo-FMS_column_picketfence

Major Update History

 - May 2021 - initial models
 - Dec 2021 - major overhaul
 - May 2022 - SW zenith angle geometric correction added
 - Jun 2022 - Bezier short char added
 - Aug 2023 - VIM added, refactored arrays
 - Nov 2026 - Ported the updated semi-grey radiation and convection schemes, split sw and lw methods, added WENO4 interpolation and n-stream DISORT benchmarks, fixed the Toon internal flux lower boundary, added band flux outputs and speed/accuracy benchmarks

Elspeth K.H. Lee - Nov 2026

This is one part of a series of codes that build upon different two-stream approaches and schemes, primarily useful for the GCM modeling community.
This is the non-grey, picket fence version, with three bands in the visible, representing incident radiation from a host star, and two bands in the IR, representing the internal radiation propagating inside the planetary atmosphere.
This is for gas giant planets with a defined internal temperature; there is no solid or liquid surface.

Some useful references for non-grey, picket fence modelling are: \
Chandrasekhar (1935) \
Parmentier & Guillot (2014) \
Parmentier et al. (2015) \
Lee et al. (2021)

The picket fence scheme builds upon the semi-grey scheme by having three visible bands and two IR bands.
The visible bands better represent the deposition of stellar flux into the atmosphere.
The IR bands represent the line and continuum opacity, with the fraction of the total flux in each band controlled by Beta_IR.
This allows cooling from two photospheric regions rather than one, giving more realistic cooling in the upper atmosphere, although isothermal profiles are still produced at very low pressures.

The fitting functions from Parmentier et al. (2014, 2015) were designed to reproduce correlated-k results and can produce realistic T-p profiles at low cost.
Small spikes can occur when the temperature-dependent opacity fits are sensitive to the profile or used outside their valid range.
Inside a GCM these are typically smoothed by the dynamical processes.

To compile enter 'make' in the src directory. To remove compiled code enter 'make clean'. \
Some compiler options for gfortran, nvfortran and intel are provided in the makefile. \
No external libraries are needed, the DISORT toolchain (including its own LAPACK routines) is included in src.

From the repository root you can also use:

```sh
make -C src -j4
./FMS_RC_picketfence
```

Run the executable from the repository root so it can read FMS_RC.nml and the hybrid sigma grids.
The supplied namelist requests one million timesteps; reduce nstep for a short test.
For nvfortran or Intel, also change DISORT_FLAGS as described in the makefile.

# Personal recommendations

For the longwave: \
'lw_DDA' - the recommended choice when scattering is present. \
'lw_Toon89' - the default in the supplied namelist. A familiar, fast scattering method. \
'lw_sc_scat_ric' - for accurate heating rates with sufficient vertical resolution.  (NOT PUBLISHED YET). \
'lw_AD' - for reference accuracy. \
'lw_AA_L' - a good simple choice if there is no longwave scattering.
'lw_TAS' - for weak scattering after checking the required accuracy.

For the shortwave: \
'sw_adding' - the default and recommended fast multiple scattering choice. \
'sw_Toon89' - the well tested Toon89 two-stream method.. \
'sw_Liou88' - for four-stream accuracy. \
'sw_AD' - for reference accuracy.  \
'sw_direct' - if there is no shortwave scattering. 

# Methods

For the shortwave method we have:

1. Direct beam only method
2. Adding method with delta-Eddington layers (Meador & Weaver 1980)
3. Spherical harmonic doubling adding method (SDA, Zhang & Li 2013)
4. Toon89 multiple scattering method
5. Two- and four-stream spherical harmonics (SH2 and SH4), converted from PICASO
6. Adding-doubling (AD) 8-stream multiple scattering, converted from AddingDoublingRT
7. Delta-four-stream discrete ordinates with adding (Liou88, Liou et al. 1988)
8. Two-stream DISORT version
9. n-stream DISORT version (benchmark)

Feautrier shortwave is in development.

For the longwave method we have:

1. Absorption approximation with exponential in tau function (AA_E)
2. Absorption approximation with linear in tau function (AA_L)
3. Short characteristics (sc) with linear interpolants and the extended absorption approximation (EAA) for scattering
4. Short characteristics 8-stream discrete ordinates with anisotropic multiple scattering, solved by Riccati invariant imbedding, with a source function technique flux integration along the Hogan (2024) angles (sc_scat_ric)
5. Toon89 multiple scattering method
6. Two- and four-stream spherical harmonics (SH2 and SH4), converted from PICASO
7. Variational iteration method (VIM), exponential and linear in tau (VIM_E and VIM_L)
8. Adding-doubling (AD) 8-stream multiple scattering, converted from AddingDoublingRT
9. Four-stream discrete ordinates adding method (DDA, Zhang et al. 2013), with an optional source function technique flux integration along the Hogan (2024) angles (use_sft, off by default)
10. Delta-two- and four-stream combination method (Fu97, Fu et al. 1997)
11. Chou et al. (1999) scaling with the Tang et al. (2018) adjustment (TAS)
12. Two-stream DISORT version
13. n-stream DISORT version (benchmark)
14. Absorption-only discontinuous finite-element method (DFE)

The direct and iterative short characteristics scattering solvers and the Feautrier method are in development and are disabled. DFE is experimental and does not include scattering.

You can also see the header comments in the source code for some additional information.

For the longwave radiation, we use the weighted essentially non oscillatory order 4 interpolation method (WENO4) from Janett et al. (2019) to interpolate the temperature from the layers to the levels (typically this is required in GCMs).
The interpolation is done in log pressure, which gives a very smooth interpolation and strong gradients in temperature are accurately interpolated through. The end points are linearly extrapolated in log(p)-log(T) space to avoid numerical issues and overshoot.

The model uses a gas giant lower boundary. Shortwave radiation is not reflected at the bottom of the column, and in the longwave the internal temperature sets the net flux at the lower boundary, F_up - F_down = sigma*Tint^4. There is no solid surface.

We also include the dry convection schemes from Exo-FMS_2: the Manabe et al. (1965) adjustment, Raymond Pierrehumbert's adjustment and mixing length theory (MLT).
These act on the column after the radiative heating of the timestep. Manabe and Pierrehumbert adjust unstable layers towards a dry adiabat, while MLT gives a time-dependent convective tendency and Kzz.

Results can be plotted using the plot_TP.py and plot_flux.py python codes in the main directory.

# Picket fence bands and boundaries

The Parmentier functions provide three visible opacity ratios (gam_V), two IR opacity ratios (gam_1 and gam_2), and the band fractions Beta_V and Beta_IR.
The chosen Freedman or Valencia Rosseland mean is multiplied by these ratios to construct each band's optical depth.
Each visible band receives Beta_V(b) * (1 - AB) * sigma*Tirr^4 of the incident stellar flux.

The complete longwave solution, including atmospheric emission and the internal flux, is multiplied by Beta_IR(b).
The lower boundary in each IR band is therefore F_up - F_down = Beta_IR(b) * sigma*Tint^4.
Summing the two IR bands gives one total internal flux, sigma*Tint^4.
The Toon port applies this condition in both the coefficient solve and the angular sweep, rather than replacing only the reported bottom net flux.

By default use_bond_fit = .True. uses the Parmentier Bond albedo fit and overrides AB, as in the previous picket fence driver.
Set use_bond_fit = .False. to use the supplied AB instead.
With explicit shortwave scattering, set use_bond_fit = .False. and AB = 0 if the solver should calculate all reflection itself.
Otherwise the prescribed Bond albedo removes a fraction of the stellar flux before the scattering calculation.

# Namelist options

In the file 'FMS_RC.nml' you can select different options that control the simulation, below are the brief descriptions and options for the parameters.

sw_scheme: \
'sw_direct' - direct beam only \
'sw_adding' - adding method \
'sw_SDA' - four-stream spherical harmonic doubling adding method \
'sw_Toon89' - Toon89 method (Exo-FMS_2 version) \
'sw_SH2' - two-stream spherical harmonic method (PICASO) \
'sw_SH4' - four-stream spherical harmonic method (PICASO) \
'sw_AD' - adding-doubling method \
'sw_Liou88' - delta-four-stream discrete ordinates method (Liou et al. 1988) \
'sw_disort_ts' - two stream disort method \
'sw_disort' - n-stream disort method (benchmark, slow)

lw_scheme: \
'lw_AA_E' - AA_E method \
'lw_AA_L' - AA_L method \
'lw_sc_linear' - sc linear method (same result as AA_L, written as short characteristics) \
'lw_sc_scat_ric' - sc 8-stream scattering method, Riccati solution \
'lw_Toon89' - Toon89 method (Exo-FMS_2 version) \
'lw_SH2' - two-stream spherical harmonic method (PICASO) \
'lw_SH4' - four-stream spherical harmonic method (PICASO) \
'lw_VIM_E' - VIM method, exponential in tau \
'lw_VIM_L' - VIM method, linear in tau \
'lw_AD' - adding-doubling method \
'lw_DDA' - four-stream discrete ordinates adding method \
'lw_Fu97' - delta-two- and four-stream combination method \
'lw_TAS' - Chou et al. (1999) scaling with the Tang et al. (2018) adjustment \
'lw_disort_ts' - two stream disort method \
'lw_disort' - n-stream disort method (benchmark, slow) \
'lw_DFE' - absorption only DFE method (requires lw_ac = 0)

opac_scheme: \
'Freedman' - Uses the Rosseland mean fitting function from Freedman et al. (2014), scaled into the five picket fence bands \
'Valencia' - Uses the Rosseland mean fitting function from Valencia et al. (2013), scaled into the five picket fence bands

The semi-grey constant, Heng and TK19 opacity modes are not included in this driver.

adj_scheme: \
'Manabe' - Manabe et al. (1965) dry convective adjustment (recommended) \
'Pierrehumbert' - Ray Pierrehumbert's dry convective adjustment ('Ray_dry' also works) \
'MLT' - Mixing length theory, also gives Kzz

The option 'none' turns off a radiation or convection scheme (e.g. to run without convective adjustment set adj_scheme = 'none'). 'None' is also accepted. An opacity scheme must still be selected.

nlay, a_sh , b_sh - the number of layers, and filenames that contain the a and b constants for the hybrid sigma grid (at least three layers for WENO interpolation)

pref - reference surface pressure (pa)

t_step - time step in seconds \
nstep - number of integer timesteps \
Rd_air - specific gas constant of the air (J kg-1 K-1)\
cp_air - specific heat capacity (constant pressure) of the air (J kg-1 K-1) \
grav - gravitational acceleration constant (m s-2) \
mu_z - cosine angle of the solar zenith angle \
Tirr - Irradiation temperature \
Tint - Internal temperature

k_V - visible opacity for the analytical initial conditions (m2 kg-1) \
k_IR - IR opacity for the analytical initial conditions (m2 kg-1) \
AB - prescribed Bond albedo, used when use_bond_fit = .False. \
use_bond_fit - use the Parmentier Bond albedo fit and override AB (.True. by default) \
fl - The Heng et al. (2011) parameter used in the analytical initial conditions \
met - metallicity in dex solar (M/H)

sw_ac(3) - shortwave single scattering albedo for each visible band (constant at all layers) \
sw_gc(3) - shortwave asymmetry factor for each visible band (constant at all layers) \
lw_ac(2) - longwave single scattering albedo for each IR band (constant at all layers) \
lw_gc(2) - longwave asymmetry factor for each IR band (constant at all layers)

iIC - Initial condition selection integer (1 = isothermal at Tint, 2 = deep Guillot approximation, 3 = Guillot (2010), 4 = Parmentier (2014/2015), 5 = Mayne (2014), 6 = adiabatic from Tdeep_init with a Tstrat_init floor) \
corr - Flag to perform the adiabatic gradient correction in the initial conditions \
table_num - Parmentier table selection (1 = with TiO and VO, 2 = without TiO and VO)

Tdeep_init - initial deep atmosphere temperature for the adiabatic initial condition \
Tstrat_init - initial temperature of stratosphere

# n-stream DISORT benchmark

sw_disort and lw_disort use the full DISORT 4.0.99 solver and are meant as accuracy benchmarks rather than for use in a GCM, they are a lot slower than the other methods.
Each call covers one band, so the benchmark cost includes three SW or two LW band solves.
The number of streams is set by the nstr parameter at the top of sw_disort_mod.f90 and lw_disort_mod.f90 (default 32), you will need to clean and recompile if you change this.

DISORT assumes a surface, which a gas giant does not have. Rather than changing DISORT, the internal flux lower boundary is found from two DISORT calls with different bottom temperatures.
Since the fluxes are linear in the bottom emission, the right combination of the two gives F_up - F_down = sigma*Tint^4 before applying the IR band weight.
The two weighted band solutions sum to the prescribed internal flux.

The DISORT toolchain is compiled in double precision (see the makefile), as single precision is not accurate enough for the small internal flux compared to the large fluxes deep in a gas giant atmosphere.

# Plotting results

plot_TP.py plots the temperature profile, heating tendencies, band optical depths and opacities.
plot_flux.py plots the broadband SW and LW fluxes and all five optical depths.
These plotting routines require NumPy and Matplotlib; plot_TP.py also imports Pandas.

The output files contain whitespace-separated numeric columns, with pressure in Pa, flux in W m-2 and temperature tendencies in K s-1.
Index 1 is the top of the column, and net flux is positive upward.

FMS_RC_ic.out - layer index, pressure and initial temperature. \
FMS_RC_pp.out - layer index, pressure, temperature, radiative tendency, convective tendency, three visible optical depths, two IR optical depths, three visible opacities, two IR opacities and Kzz. The first 15 columns retain their previous order. \
FMS_RC_flx.out - level index, pressure, broadband SW up/down/net, broadband LW up/down/net, three visible optical depths and two IR optical depths. \
FMS_RC_band_flx.out - level index, pressure, SW up (3), SW down (3), SW net (3), LW up (2), LW down (2) and LW net (2). These fluxes already include the band weights.

The fluxes and tendencies describe the final radiation step, before its temperature update, as in the semi-grey driver.
Output files are replaced on each run.

# Gaussian Ordinance values

In methods that use Gaussian quadrature to perform the mu integration in intensity (to calculate the flux), various values and weights can be changed for testing at the top of the modules.
Most longwave methods use the 4 node Gauss-Jacobi quadrature from Hogan (2024).
The schemes where the quadrature points are also the scattering streams (sc_scat_ric, AD) use Gauss-Legendre points. sc_scat_ric uses its Gauss-Legendre streams for the scattering only, the fluxes are integrated along the Hogan (2024) angles (use_sft, on by default).
You will need to clean and recompile the code if these are changed. Beware, some codes use the classical Gauss-Legendre method (flux = 2pi * intensity), while others include the weighting in the intensity (flux = pi * intensity).

# Future developments

Additional opacity schemes from the literature can be added. \
Window fractions and functions.
