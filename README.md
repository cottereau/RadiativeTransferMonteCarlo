# Radiative Transfer Monte Carlo

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Monte Carlo approximation of the radiative transfer equation for acoustic and elastic waves in random media. This code computes wave energy evolution and propagation in complex scattering environments, with applications to seismology, ultrasonics, and wave physics.

## What It Does

This code solves radiative transfer equations using Monte Carlo particle tracing. It is designed for acoustic and elastic wave-energy transport in random scattering media.
- **Acoustic and elastic wave propagation** in both two- and three-dimensional media.
- **Bounded and unbounded media** with multiple scattering.
- **Wave energy density evolution** as a function of space, propagation direction, and time.
- **Validated against published results** in the radiative transfer literature.

The code supports multiple geometries (Cartesian, spherical, cylindrical) and wave types (acoustic, elastic), making it suitable for seismic wave modeling, ultrasonic testing, and wave physics research.

## Why Use This Code

- **Validated approach**: Results compared against five key papers in radiative transfer theory
- **Flexible geometry**: 2D/3D isotropic and anisotropic scattering operators
- **Efficient computation**: Parallel processing support via MATLAB Parallel Computing Toolbox
- **Research-ready**: Used in peer-reviewed publications and academic research

## Quick Start

### Prerequisites

- MATLAB R2016b or later
- (Optional) MATLAB Parallel Computing Toolbox for faster execution

### Installation

```matlab
% Option 1: Add to path with full path argument
AddRTMCLib('/path/to/RadiativeTransferMonteCarlo')

% Option 2: Run from within the repository folder
cd RadiativeTransferMonteCarlo
AddRTMCLib
```

The `AddRTMCLib.m` script adds the repository to your MATLAB path, allowing access to functions from any working directory.

### Basic Example: Acoustic Wave Propagation

See [Examples/mainAcoustics.m](Examples/mainAcoustics.m) for a complete example. Here is the essential workflow:

```matlab
% Note : The following example uses a simple consistent unit system for illustration 
% in physical simulations, all lengths, velocities, times, frequencies, and 
% correlation lengths must be expressed in compatible units.

% Define geometry: 3D unbounded domain, spherical observation frame
geometry = struct('dimension', 3, 'frame', 'spherical');

% Define initial energy distribution: isotropic point source at the origin
source = struct('numberParticles', 1e6, ...
                'type', 'point', ...
                'position', [0 0 0], ...
                'direction', 'uniform', ...
                'lambda', 0.05);

% Define observation grid and time
% In the spherical frame: x = radius, y = azimuth, z = elevation.
% We resolve the radial coordinate and integrate over angles and directions.
observation = struct('x', 0:0.05:6, ...
                     'y', [-pi pi], ...
                     'z', [-pi/2 pi/2], ...
                     'directions', [0 pi], ...
                     'time', 0:0.02:5);

% Create acoustic material
freq = 10;          % frequency
acoustics = true;   % acoustic waves
v = 1;              % sound speed
cv = [0.05 0.05];   % coefficients of variation of compressibility and density
corr = 0.0;         % correlation coefficient between compressibility and density
acf = 'exp';        % autocorrelation model
lc = 0.1;           % correlation length

material = MaterialClass(geometry, freq, acoustics, v, cv, corr, acf, lc);

% No intrinsic attenuation
material.Q = Inf; % default value

% Run simulation
obs = radiativeTransfer(geometry, source, material, observation);

% Check the total recorded energy
figure;
plot(obs.t, obs.energy, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Recorded energy');
title('Energy conservation check');
grid on;

% Check the energy density
figure;
semilogy(obs.t, squeeze(obs.energyDensity(10,:,:)), 'LineWidth', 1.5);
xlabel('Time');
ylabel('Recorded energy density');
title('Energy density at a given point');
grid on;
```

## Project Structure

```
├── Examples/                 # Example scripts and validation cases
│   ├── mainAcoustics.m       # Acoustic wave examples
│   ├── mainElastics.m        # Elastic wave examples
│   └── RunAllTests.m         # Run all validation comparisons
├── @MaterialClass/           # Material properties and scattering cross-sections
├── +Comparison/              # Comparison functions against literature
├── AddRTMCLib.m              # Setup and path configuration
└── radiativeTransfer.m       # Main solver
```

## Main Functions

| Function | Purpose |
|----------|---------|
| `radiativeTransfer()` | Solver for bounded and unbounded domains |
| `MaterialClass()` | Material specification and scattering properties |
| `Examples/RunAllTests` | Run validation against published results |

## Input Data Format

The main solver is called as

```matlab
obs = radiativeTransfer(geometry, source, material, observation);
```

where `geometry`, `source`, and `observation` are MATLAB structures, and `material` is a `MaterialClass` object.

### Geometry

```matlab
geometry = struct('dimension', 3, 'frame', 'cylindrical');
```

Available fields:

| Field       | Type                   | Description                                                                                                    |
| ----------- | ---------------------- | -------------------------------------------------------------------------------------------------------------- |
| `dimension` | integer                | Spatial dimension, `2` or `3`. For elastic 2D simulations, the implemented model corresponds to the in-plane P-SV case.|
| `frame`     | char                   | Observation frame: `'spherical'` (default), `'cylindrical'`, or `'cartesian'`. |
| `bnd`       | struct array (optional)| Boundary definitions. By default the domain is unbounded.                                                     |

The boundary structure has the form

```matlab
geometry.bnd(i) = struct('dir', dir, 'val', val, 'type', type);
```

where:

| Field     | Meaning                                                                     |
| --------- | --------------------------------------------------------------------------- |
| `dir = 1` | Plane boundary normal to the `x` direction, located at `x = val`.           |
| `dir = 2` | Plane boundary normal to the `y` direction, located at `y = val`.           |
| `dir = 3` | Plane boundary normal to the `z` direction, located at `z = val`.           |
| `dir = 4` | Cylindrical boundary of radius `val` around the `z` axis.                   |
| `val`     | Boundary coordinate or radius.                                              |
| `type`    | `'reflective'` (default) or `'absorbing'`. |

Example:

```matlab
geometry.bnd(1) = struct('dir', 3, 'val', 0);
```

defines a reflective plane boundary at `z = 0`.

### Source

Available fields:

| Field             | Type                            | Description                                                                                                 |
| ----------------- | ------------------------------- | ----------------------------------------------------------------------------------------------------------- |
| `numberParticles` | integer                         | Number of Monte Carlo particles.                                                                            |
| `type`            | char                            | `'point'` (default) or `'plane'`                                                                 |
| `position`        | vector                          | Source position in Cartesian coordinates. In 2D, `[x y]` is accepted and internally converted to `[x y 0]`. |
| `lambda`          | float                           | Width of the Gaussian initial spatial distribution.                                                         |
| `direction` | char, numeric vector, or scalar axis index | For point sources, use `'uniform'` (default), `'outgoing'`, or `'upper'`. For plane sources, use a 3-component vector such as `[0 0 1]`. Axis indices `±1`, `±2`, `±3` are equivalent to the corresponding Cartesian unit vectors. |
| `radial`          | function/distribution, optional | User-defined radial distribution for point sources.                                                         |
| `extent` | scalar or vector | Plane-source aperture size. For a rectangular aperture, use `[L1 L2]`, where `L1` and `L2` are the side lengths in the plane perpendicular to `source.direction`. For a circular aperture, use a scalar radius `R` with `source.aperture = 'circle'`. |
| `aperture`        | char, optional                  | Plane-source aperture: `'rectangle'` (default) or `'circle'`.                                 |
| `polarization` | char, elastic only | Initial wave mode for elastic simulations. Use `'P'` (default) for compressional-wave particles or `'S'` for shear-wave particles. |

For a point source, `source.direction` can be:

| Value                  | Meaning                                                                                                    |
| ---------------------- | ---------------------------------------------------------------------------------------------------------- |
| `'uniform'` (default) | Isotropic source direction.                                                                                |
| `'outgoing'`           | Initial propagation direction is aligned with the initial particle position relative to the source center. |
| `'upper'`              | Isotropic source over the upper hemisphere.                                                                |

### Observation

The observation structure defines the histogram bins used to record particle energy.

| Field        | Description                                                                               |
| ------------ | ----------------------------------------------------------------------------------------- |
| `x`          | First spatial bin coordinate. Its meaning depends on `geometry.frame`.                    |
| `y`          | Second spatial/angular bin coordinate.                                                    |
| `z`          | Third spatial/angular bin coordinate, used in 3D.                                  |
| `directions` | Bins for the angle between the propagation direction and the position vector. |
| `time`       | Observation times.                                                                        |

Only two variables among `x`, `y`, `z`, and `directions` should have more than one bin, because the code stores two-dimensional histograms. `directions` is an angle in `[0, pi]`. Use `directions = [0 pi]` when you want to integrate over all propagation directions.

Frame conventions:

| Frame           | Meaning of `x` | Meaning of `y`               | Meaning of `z`                     |
| --------------- | -------------- | ---------------------------- | ---------------------------------- |
| `'cartesian'`   | Cartesian `x`  | Cartesian `y`                | Cartesian `z`                      |
| `'cylindrical'` | Radius `r`     | Azimuth angle in `[-pi, pi]` | Cartesian `z`                      |
| `'spherical'`   | Radius `r`     | Azimuth angle in `[-pi, pi]` | Elevation angle in `[-pi/2, pi/2]` |

### Material

Materials are created with

```matlab
material = MaterialClass(geometry, freq, acoustics, v, cv, corr, acf, lc);
```

where:

| Argument / field | Description                                                                                                                            |
| ---------------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| `geometry`       | Geometry structure.                                                                                                                    |
| `freq`           | Frequency in Hz.                                                                                                                       |
| `acoustics`      | `true` for acoustic waves, `false` for elastic waves.                                                                                  |
| `v`              | Sound speed for acoustics, or `[vp vs]` for elastics.                                                                                  |
| `cv`             | Coefficients of variation of the fluctuating material parameters.                                                                      |
| `corr`           | Correlation coefficient(s) between fluctuating material parameters.                                                                    |
| `acf`            | Autocorrelation model.                                                                                                                 |
| `lc`             | Correlation length.                                                                                                                    |
| `material.Q` | Intrinsic quality factor. Default: `Inf`, meaning no intrinsic attenuation. For acoustics, use a scalar. For elastic waves, use either a scalar value for both modes or `[QP QS]`. |

For acoustic waves,

```matlab
cv   = [CV_kappa, CV_rho]; % coefficients of variation of compressibility and density
corr = corr_kappa_rho;     % correlation coefficient between compressibility and density
```

where `kappa` and `rho` are compressibility and mass density, respectively.

For elastic waves,

```matlab
cv   = [CV_lambda, CV_mu, CV_rho];                     % coefficients of variation of Lamé parameters and density
corr = [corr_lambda_mu, corr_lambda_rho, corr_mu_rho]; % mutual correlation coefficients
```

where `lambda` and `mu` are Lamé coefficients. The code uses direct fluctuations of `lambda`, `mu`, and `rho`.

### Autocorrelation and Power Spectral Density Models

The argument `acf` selects the autocorrelation / power spectral density model used to compute the differential scattering cross-sections.

Currently supported options include:

| `acf` value            | Description                                                         |
| ---------------------- | ------------------------------------------------------------------- |
| `'exp'`                | Exponential correlation model.                                      |
| `'gaussian'`           | Gaussian correlation model.                                         |
| `'power_law'`          | Power-law correlation model.                                        |
| `'VonKarman'`          | von Kármán correlation model. Requires `material.SpectralParam.nu`. |
| `'triangular'`         | Triangular spectral model.                                          |
| `'low_pass'`           | Low-pass spectral model.                                            |
| `'monodispersesphere'` | Microstructure-inspired model based on mono-disperse inclusions.    |
| `'image'`              | Spectrum estimated from an image or binary microstructure.          |
| `'Imported'`           | User-imported autocorrelation function.                             |

### Autocorrelation Normalization

The code assumes that the random medium is statistically homogeneous and isotropic. Therefore, the normalized autocorrelation function depends only on the lag distance `r`:

```text
R = R(r),     with R(0) = 1.
```

The analytical ACF models are written in terms of the normalized lag distance `s = r / Lc` where `Lc` is the correlation length. The correlation length is defined from the dimensional ACF as follows:

```text
1D:  Lc   = 2 * int_0^inf R(r) dr

2D:  Lc^2 = 2 * int_0^inf r * R(r) dr

3D:  Lc^3 = 3 * int_0^inf r^2 * R(r) dr
```

Equivalently, if the ACF is written as a dimensionless function of `s = r/Lc`, then the normalized models satisfy

```text
1D:  1 = 2 * int_0^inf R(s) ds

2D:  1 = 2 * int_0^inf s * R(s) ds

3D:  1 = 3 * int_0^inf s^2 * R(s) ds
```

### Intrinsic Attenuation

Intrinsic attenuation can be included through the quality factor `material.Q`.

```matlab
material.Q = Inf;       % default: no intrinsic attenuation
material.Q = 100;       % for acoustic waves
material.Q = [200 100]; % for elastic waves: [QP QS]
```

Intrinsic attenuation is modeled by stochastic killing of particles during propagation. For a time step `dt`, particles are removed with probability

```text
pKill = 1 - exp(-omega*dt/Q)
```

where `omega = 2*pi*material.Frequency` is the angular frequency.

## Theoretical Background

The code implements the radiative transfer equation approach described in:

- **Ryzhik, L., Papanicolaou, G. and Keller, J.B.**, 1996. *Transport equations for elastic and other waves in random media*. Wave Motion, 24(4), pp. 327-370. [DOI](https://doi.org/10.1016/S0165-2125(96)00021-2)

Results are validated against these key papers:

1. **Paasschens, J.C.J.**, 1997. *Solution of the time-dependent Boltzmann equation*. Physical Review E, 56(1), pp. 1135-1141. [DOI](https://doi.org/10.1103/PhysRevE.56.1135)
2. **Yoshimoto, K.**, 2000. *Monte Carlo simulation of seismogram envelopes in scattering media*. Journal of Geophysical Research: Solid Earth, 105(B3), pp. 6153-6161. [DOI](https://doi.org/10.1029/1999JB900437)
3. **Nakahara, H. and Yoshimoto, K.**, 2011. *Radiative transfer of elastic waves in two-dimensional isotropic scattering media: Semi-analytical approach for isotropic source radiation*. Earth, Planets and Space, 63, pp. 459-468. [DOI](https://doi.org/10.5047/eps.2011.03.006)
4. **Hoshiba, M.**, 1991. *Simulation of multiple-scattered coda wave excitation based on the energy conservation law*. Physics of the Earth and Planetary Interiors, 67(1-2), pp. 123-136. [DOI](https://doi.org/10.1016/0031-9201(91)90066-Q)
5. **Sato, H.**, 1994. *Multiple isotropic scattering model including P-S conversions for the seismogram envelope formation*. Geophysical Journal International, 117(2), pp. 487-494. [DOI](https://doi.org/10.1111/j.1365-246X.1994.tb03946.x)

### Monte Carlo Method References

The implementation is based on the Monte Carlo method. Additional details on this methodology can be found in:

- **Gomez, C. and Pinaud, O.**, 2018. *Monte Carlo methods for radiative transfer with singular kernels*. SIAM Journal on Scientific Computing, 40(3), pp. A1714-A1741. [DOI](https://doi.org/10.1137/17M1134755)
- **Lapeyre, B., Pardoux, E. and Sentis, R.**, 2003. *Introduction to Monte Carlo methods for transport and diffusion equations*. Vol. 6. OUP Oxford. [DOI](https://doi.org/10.1093/oso/9780198525929.001.0001)
- **Margerin, L., Campillo, M. and Van Tiggelen, B.**, 2000. *Monte Carlo simulation of multiple scattering of elastic waves*. Journal of Geophysical Research: Solid Earth, 105(B4), pp. 7873-7892. [DOI](https://doi.org/10.1029/1999JB900359)

## Publications

- **Corrêa, L.D.A., Khazaie, S., Gomez, C. and Cottereau, R.**, 2026. *Quantitative error assessment of radiative transfer approximations of acoustic wave energies in unbounded and bounded random media*. Wave Motion, p.103715. [DOI](https://doi.org/10.1016/j.wavemoti.2026.103715)

## How to Cite

If you use this code in your research, please cite the following paper:

```bibtex
@article{CorreaKhazaieGomezCottereau2026,
  title   = {Quantitative error assessment of radiative transfer approximations of acoustic wave energies in unbounded and bounded random media},
  author  = {Corr{\^e}a, L. D. A. and Khazaie, S. and Gomez, C. and Cottereau, R.},
  journal = {Wave Motion},
  pages   = {103715},
  year    = {2026},
  doi     = {10.1016/j.wavemoti.2026.103715}
}
```

This helps us track the use of the code and supports further development.

## Support & Documentation

- **Examples**: See [Examples/](Examples/) directory for complete workflows
- **Validation**: Run [Examples/RunAllTests.m](Examples/RunAllTests.m) to reproduce published results

## Contributing & Support

For questions, bug reports, or contributions:

- **Lead Maintainer:** [Régis Cottereau](mailto:cottereau@lma.cnrs-mrs.fr) — Laboratoire de Mécanique et d'Acoustique (LMA), CNRS, Marseille
- **Contributors:** 
  - [Lúcio de Abreu Corrêa](mailto:de-abreu-correa@lma.cnrs-mrs.fr) — LMA, CNRS
  - [Shahram Khazaie](mailto:Shahram.Khazaie@univ-nantes.fr) — Institut de Recherche en Génie Civil et Mécanique (GeM), Nantes Université

## License

This project is licensed under the **GNU General Public License v3.0** — see [LICENSE](LICENSE) for details.
