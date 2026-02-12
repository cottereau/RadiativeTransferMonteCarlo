# Radiative Transfer Monte Carlo

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Monte Carlo approximation of the radiative transfer equation for acoustic and elastic waves in random media. This project computes wave energy evolution and propagation in complex scattering environments, with applications to seismology, ultrasonics, and wave physics.

## What It Does

This project solves the radiative transfer equation using Monte Carlo particle tracing methods. It's designed for:

- **Acoustic and elastic wave propagation** in 2D and 3D domains
- **Bounded and unbounded media** with random scattering
- **Wave energy density evolution** across space, direction, and time
- **Validated against published results** in the radiative transfer literature

The code supports multiple geometries (Cartesian, spherical, cylindrical) and wave types, making it suitable for seismic wave modeling, ultrasonic testing, and wave physics research.

## Why Use This Project

- **Validated approach**: Results compared against five key papers in radiative transfer theory
- **Flexible geometry**: 2D/3D isotropic and anisotropic scatter operators
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

See [Examples/mainAcoustics.m](Examples/mainAcoustics.m) for a complete example. Here's the essential workflow:

```matlab
% Define geometry (3D unbounded domain)
geometry = struct('dimension', 3, 'frame', 'cylindrical');

% Define initial energy distribution (plane source)
source = struct('numberParticles', 1e6, ...
                'type', 'plane', ...
                'position', [0 0 -2], ...
                'direction', 3, ...
                'extent', [10 10], ...
                'lambda', 0.1);

% Define observation grid and time
observation = struct('x', -2:.1:2, ...
                     'y', linspace(-pi,pi,30), ...
                     'z', [-Inf Inf], ...
                     'directions', [-Inf Inf], ...
                     'time', 0:.1:5);

% Create material (isotropic acoustic medium)
material = MaterialClass(3, 'isotropic', 'acoustic');
material.rho = 2000;          % density (kg/m³)
material.v = 2000;            % sound speed (m/s)
material.Frequency = 100;     % frequency (Hz)
material.correlation_coefficients = [0.3 0.5 0.3];

% Run simulation
obs = radiativeTransferUnbounded(geometry, source, material, observation);
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

The solver accepts four structured arrays:

### Geometry

- `dimension` (int): 2 or 3
- `frame` (char, optional): `'spherical'` (default), `'cartesian'`, or `'cylindrical'`
- `bnd` (struct array, optional): Boundary definitions for bounded domains

### Source

- `numberParticles` (int): Number of Monte Carlo particles
- `type` (char): `'point'` or `'plane'`
- `position` (1×3 float): Source location (always Cartesian)
- `lambda` (float): Gaussian source width (controls wavelength)
- `direction` (point source: char `'uniform'`/`'outgoing'`; plane source: int ±1/±2/±3)
- `extent` (plane source only): Source dimensions [Δy, Δz]

### Observation

- `x`, `y`, `z`, `directions` (N×1 vectors): Spatial and angular bins
- `time` (M×1 vector): Observation times

For spherical coordinates (default): x=radius, y=azimuth [−π, π], z=elevation [−π/2, π/2]

### Material

Created via `MaterialClass()`. Defines medium type (isotropic/anisotropic), wave type (acoustic/elastic), and physical properties:

```matlab
material = MaterialClass(3, 'isotropic', 'acoustic');
material.rho = 2000;                    % density
material.v = 2000;                      % sound speed
material.Frequency = 100;               % frequency
material.correlation_coefficients = [0]; % scattering correlation
```

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

## Support & Documentation

- **API Documentation**: Detailed input/output specifications are embedded in function headers.
- **Examples**: See [Examples/](Examples/) directory for complete workflows
- **Validation**: Run [Examples/RunAllTests.m](Examples/RunAllTests.m) to reproduce published results

## Contributing & Support

For questions, bug reports, or contributions:

- **Lead Maintainer:** [Régis Cottereau](mailto:cottereau@lma.cnrs-mrs.fr) — Laboratoire de Mécanique et d'Acoustique (LMA), CNRS, Marseille
- **Contributors:** 
  - [Lucio de Abreu Corrêa](mailto:de-abreu-correa@lma.cnrs-mrs.fr) — LMA, CNRS
  - [Shahram Khazaie](mailto:Shahram.Khazaie@univ-nantes.fr) — Institut GeM, Nantes University

## License

This project is licensed under the **GNU General Public License v3.0** — see [LICENSE](LICENSE) for details.
