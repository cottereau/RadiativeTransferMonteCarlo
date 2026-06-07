# Radiative Transfer Monte Carlo

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Monte Carlo approximation of the radiative transfer equation for acoustic and elastic waves in random media. This code computes wave energy evolution and propagation in complex scattering environments, with applications to seismology, ultrasonics, and wave physics.

## What It Does

This code solves radiative transfer equations using Monte Carlo particle tracing. It is designed for acoustic and elastic wave-energy transport in random scattering media.
- **Acoustic and elastic wave propagation** in 2D and 3D domains.
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
v = 1;              % sound speed, dimensionless units
cv = [0.05 0.05];   % [CV_kappa, CV_rho]
corr = 0.0;         % corr(kappa,rho)
acf = 'exp';        % autocorrelation model
lc = 0.1;           % correlation length

material = MaterialClass(geometry, freq, acoustics, v, cv, corr, acf, lc);

% No intrinsic attenuation
material.Q = Inf;

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
├── Examples/                    # Example scripts and validation cases
│   ├── mainAcoustics.m         # Acoustic wave examples
│   ├── mainElastics.m          # Elastic wave examples
│   └── RunAllTests.m           # Run all validation comparisons
├── @MaterialClass/             # Material properties and scattering cross-sections
├── +Comparison/                # Comparison functions against literature
├── AddRTMCLib.m               # Setup and path configuration
└── radiativeTransfer.m        # Main solver
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
| `dimension` | integer                | Spatial dimension, `2` or `3`.                                                                                 |
| `frame`     | char                   | Observation frame: `'spherical'`, `'cylindrical'`, or `'cartesian'`. The default is `'spherical'`. |
| `bnd`       | struct array, optional | Boundary definitions. If omitted, the domain is unbounded.                                                     |

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
| `type`    | `'reflective'` or `'absorbing'`. The default is `'reflective'`. |

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
| `type`            | char                            | `'point'` or `'plane'`. Default: `'point'`.                                                                 |
| `position`        | vector                          | Source position in Cartesian coordinates. In 2D, `[x y]` is accepted and internally converted to `[x y 0]`. |
| `lambda`          | float                           | Width of the Gaussian initial spatial distribution.                                                         |
| `direction`       | char, integer, or vector        | Direction option. Its meaning depends on `source.type`.                                                     |
| `radial`          | function/distribution, optional | User-defined radial distribution for point sources.                                                         |
| `extent`          | scalar or vector                | Plane-source aperture size.                                                                                 |
| `aperture`        | char, optional                  | Plane-source aperture: `'rectangle'` or `'circle'`. Default: `'rectangle'`.                                 |
| `polarization`    | char, elastic only              | `'P'` or `'S'`. Default: `'P'`.                                                                             |

For a point source, `source.direction` can be:

| Value                  | Meaning                                                                                                    |
| ---------------------- | ---------------------------------------------------------------------------------------------------------- |
| omitted or `'uniform'` | Isotropic source direction.                                                                                |
| `'outgoing'`           | Initial propagation direction is aligned with the initial particle position relative to the source center. |
| `'upper'`              | Isotropic source over the upper hemisphere.                                                                |

For a plane source, `source.direction` can be either a vector or an axis index:

| Value              | Meaning                         |
| ------------------ | ------------------------------- |
| `[1 0 0]` or `1`   | Propagation along positive `x`. |
| `[-1 0 0]` or `-1` | Propagation along negative `x`. |
| `[0 1 0]` or `2`   | Propagation along positive `y`. |
| `[0 -1 0]` or `-2` | Propagation along negative `y`. |
| `[0 0 1]` or `3`   | Propagation along positive `z`. |
| `[0 0 -1]` or `-3` | Propagation along negative `z`. |

For rectangular plane sources, `extent = [L1 L2]` gives the two transverse aperture lengths. For circular plane sources, use

```matlab
source.aperture = 'circle';
source.extent = R;
```

where `R` is the disk radius.

### Observation

The observation structure defines the histogram bins used to record particle energy.

| Field        | Description                                                                               |
| ------------ | ----------------------------------------------------------------------------------------- |
| `x`          | First spatial bin coordinate. Its meaning depends on `geometry.frame`.                    |
| `y`          | Second spatial/angular bin coordinate.                                                    |
| `z`          | Third spatial/angular bin coordinate, used mainly in 3D.                                  |
| `directions` | Bins for the angle between the propagation direction and the observation position vector. |
| `time`       | Observation times.                                                                        |

Only two variables among `x`, `y`, `z`, and `directions` should have more than one bin, because the code stores two-dimensional histograms. `directions` is an angle in `[0, pi]`. Use `directions = [0 pi]` when you want to integrate over all propagation directions.

Frame conventions:

| Frame           | Meaning of `x` | Meaning of `y`               | Meaning of `z`                     |
| --------------- | -------------- | ---------------------------- | ---------------------------------- |
| `'cartesian'`   | Cartesian `x`  | Cartesian `y`                | Cartesian `z`                      |
| `'cylindrical'` | Radius `r`     | Azimuth angle in `[-pi, pi]` | Cartesian `z`                      |
| `'spherical'`   | Radius `r`     | Azimuth angle in `[-pi, pi]` | Elevation angle in `[-pi/2, pi/2]` |

In 3D, `directions` is an angle in `[0, pi]`. Internally, the code uses the corresponding `cos(psi)` bins.

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

By default,

```matlab
material.Q = Inf;
```

so no intrinsic attenuation is applied.

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

The code uses dimension-dependent normalized spectra. In `d` dimensions, the dimensional spectrum scales as

```text
Rhat_d(q) = variance * Lc^d * Phi_d(q*Lc)
```

Therefore, 2D simulations use an (Lc^2) scaling, while 3D simulations use an (Lc^3) scaling.

For elastic 2D simulations, the implemented built-in model corresponds to the in-plane P-SV case. The 2D elastic formulas use the 2D power spectral density and a single in-plane SV polarization.

### Intrinsic Attenuation

Intrinsic attenuation can be included through the quality factor `material.Q`.

```matlab
material.Q = Inf;       % default: no intrinsic attenuation
material.Q = 100;       % same Q for all waves
material.Q = [200 100]; % elastic case: [QP QS]
```

Intrinsic attenuation is modeled by stochastic killing of Monte Carlo particles during propagation. For a time step `dt`, particles are removed with probability

```text
pKill = 1 - exp(-omega*dt/Q)
```

where `omega = 2*pi*material.Frequency` is the angular frequency.

For acoustic simulations, use a scalar `Q`. For elastic simulations, either a scalar value can be used for both P and S waves, or a two-component vector `[QP QS]` can be used.

## Theoretical Background

The code implements the radiative transfer equation approach described in:

- **Ryzhik, Papanicolaou & Keller (1996)**: Foundational theory for elastic/acoustic waves in random media. [DOI](https://doi.org/10.1016/S0165-2125(96)00021-2)

Results are validated against these key papers:

1. J. C. J. Paasschens. Solution of the time-dependent Boltzmann equation, *Phys. Rev. E* 56(1), pp. 1135-1141 (1997). [DOI](https://doi.org/10.1103/PhysRevE.56.1135)
2. M. Hoshiba. Simulation of multiple-scattered coda wave excitation based on the energy conservation law. *Phys. Earth Planet. Int.* 67, pp. 123-136 (1991). [DOI](https://doi.org/10.1016/0031-9201(91)90066-Q)
3. H. Nakahara, K. Yoshimoto. Radiative transfer of elastic waves in two-dimensional isotropic scattering media: semi-analytical approach for isotropic source radiation. *Earth Planets Space* 63, pp. 459-468 (2011). [DOI](https://doi.org/10.5047/eps.2011.03.006)
4. H. Sato. Multiple isotropic scattering model including P-S conversions for the seismogram envelope formation. *Geophys. J. Int* 117, pp. 487-494 (1994). [DOI](https://doi.org/10.1111/j.1365-246X.1994.tb03946.x)
5. K. Yoshimoto, Monte Carlo simulation of seismogram envelopes in scattering media. *J. Geophys. Res.: Solid Earth* (2000). [DOI](https://doi.org/10.1029/1999JB900437)

### Monte Carlo Method References

The implementation is based on the Monte Carlo method. Additional details on this methodology can be found in:

- **C. Gomez & O. Pinaud (2018).** Monte Carlo methods for radiative transfer with singular kernels. *SIAM J. Sci. Comp.*, 40(3), pp. A1714-A1741. [DOI](https://doi.org/10.1137/17M1134755)
- **B. Lapeyre, E. Pardoux & R. Sentis.** Introduction to Monte-Carlo Methods for Transport and Diffusion Equations.
- **L. Margerin, M. Campillo & B.A. van Tiggelen (2000).** Monte Carlo simulation of multiple scattering of elastic waves. *J. Geophys. Res.*, 105(B4), pp. 7873-7892. [DOI](https://doi.org/10.1029/1999JB900296)

## Publications

- **Corrêa, L. D. A., Khazaie, S., Gomez, C., & Cottereau, R. (2026).** Quantitative error assessment of radiative transfer approximations of acoustic wave energies in unbounded and bounded random media. *Wave Motion*, 103715. [DOI](https://doi.org/10.1016/j.wavemoti.2026.103715)

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

- **API Documentation**: Detailed input/output specifications are embedded in function headers.
- **Examples**: See [Examples/](Examples/) directory for complete workflows
- **Validation**: Run [Examples/RunAllTests.m](Examples/RunAllTests.m) to reproduce published results

## Contributing & Support

For questions, bug reports, or contributions:

- **Lead Maintainer:** [Régis Cottereau](mailto:cottereau@lma.cnrs-mrs.fr) — Laboratoire de Mécanique et d'Acoustique (LMA), CNRS, Marseille
- **Contributors:** 
  - [Lucio de Abreu Corrêa](mailto:de-abreu-correa@lma.cnrs-mrs.fr) — LMA, CNRS
  - [Shahram Khazaie](mailto:Shahram.Khazaie@univ-nantes.fr) — Institut de Recherche en Génie Civil et Mécanique (GeM), Nantes Université

## License

This project is licensed under the **GNU General Public License v3.0** — see [LICENSE](LICENSE) for details.
