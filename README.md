# Radiative Transfer Monte Carlo

## Introduction

This project proposes an approximation of solutions to the radiative transfer equation using the Monte Carlo method. It provides solutions for acoustic waves and elastic waves, both in 2D and 3D, and in bounded and unbounded media. The code is written in MATLAB. It is primarily developed at the Laboratoire de Mécanique et d'Acoustique (LMA) in Marseille, France and Institut de Recherche en Génie Civil et Mécanique (GeM), Nantes University Nantes France.

- **Contact:** [Régis Cottereau](mailto:cottereau@lma.cnrs-mrs.fr)
- **Contributors:**  [Régis Cottereau](mailto:cottereau@lma.cnrs-mrs.fr),  [Lucio de Abreu Corrêa](mailto:de-abreu-correa@lma.cnrs-mrs.fr), [Shahram Khazaie](mailto:Shahram.Khazaie@univ-nantes.fr)

## Installation

The installation can be done using the `AddRTMCLib.m` script. This function will add the folder containing the Radiative Transfer Monte Carlo repository to the MATLAB path, allowing access to the repository functions from any work folder. 

There are two ways to add the library:

1. Run `AddRTMCLib(baseFolder)`, where `baseFolder` is the folder path of the repository.
2. Run `AddRTMCLib` directly from within the repository folder.

## Literature Comparison and Validation
The theory of radiative transfer is (in our opinion) well described in the paper:
- Ryzhik, L., Papanicolaou, G., & Keller, J. B. (1996). Transport equations for elastic and other waves in random media. *Wave Motion*, 24(4), 327-370. [DOI](https://doi.org/10.1016/S0165-2125(96)00021-2)

The approximations provided by this code can be compared to the results described in the following papers:
1. J. C. J. Paasschens. Solution of the time-dependent Boltzmann equation, *Phys. Rev. E* 56(1), pp. 1135-1141 (1997). [DOI](https://doi.org/10.1103/PhysRevE.56.1135)
2. M. Hoshiba. Simulation of multiple-scattered coda wave excitation based on the energy conservation law. *Phys. Earth Planet. Int.* 67, pp. 123-136 (1991). [DOI](https://doi.org/10.1016/0031-9201(91)90066-Q)
3. H. Nakahara, K. Yoshimoto. Radiative transfer of elastic waves in two-dimensional isotropic scattering media: semi-analytical approach for isotropic source radiation. *Earth Planets Space* 63, pp. 459-468 (2011). [DOI](https://doi.org/10.5047/eps.2011.03.006)
4. H. Sato. Multiple isotropic scattering model including P-S conversions for the seismogram envelope formation. *Geophys. J. Int* 117, pp. 487-494 (1994). [DOI](https://doi.org/10.1111/j.1365-246X.1994.tb03946.x)
5. K. Yoshimoto, Monte Carlo simulation of seismogram envelopes in scattering media. *J. Geophys. Res.: Solid Earth* (2000)

The function `mainLiteratureComparison(type)` can be used to run the comparison between our code and the literature presented above. The input of the function is a chain of characters in the following list: 'all', '2dIsotropicAcoustic', '3dIsotropicAcoustic', '2dAnisotropicAcoustic', '3dAnisotropicAcoustic', '2dIsotropicElastic', '3dIsotropicElastic', or '3dAnisotropicElastic' .

## Usage

Two main examples are provided in this repository: `mainAcoustics.m` and `mainElastics.m`. These examples are tailored for solving the acoustic and elastic wave equations, respectively.

The main routine is `radiativeTransferUnbounded`. It computes an approximation of the radiative transfer solution using the specified input parameters and returns the results in a structured array (named `obs` below).

Usage: `obs = radiativeTransferUnbounded( geometry, source, material, observation );`

### Input Data

Here is the format for the four structured arrays used as input: `source`, `material`, `observation`, `geometry`:

#### Geometry
  - `dimension` (integer 2 or 3): dimensionality of the propagation space.
  - `frame` (chain of characters 'spherical' (default) or 'cartesian'): frame in which coordinates are expressed on output.
  - `bnd` (structured array(s)): description of the boundaries through their normal (field `dir`, with values 1=`x`, 2=`y`, 3=`z`) and position (indicated by field `pos`).

#### Source
Source is actually an initial distribution of the energy in phase-space. The shape is controlled by the fields below:
- `numberParticles` (integer): number of particles used to discretized the initial distribution of the energy
- `type` (chain of characters 'plane' or 'point'): energy is initially located around a plane or a point
- `position` (1x3 float vector): central position around which the initial energy is distributed
- `lambda` (scalar float): standard deviation of the Gaussian source. This controls the wavelength/wavenumber simulated, and is related to frequency.
- `direction` (with point source, chain of characters 'outgoing' or 'uniform' ): initial energy radiates away from the center position or randomly
- `radial` (with point source, function handle): describes the initial distribution in space (and overrules 'lambda')
- `direction` (with plane source, signed integer): initial direction of the energy 1=+x, 2=+y, 3=+z, -1=-x, -2=-y, -3=-z 
- `extent` (with plane source, 1x2 float vector): dimensions of the plane perpendicular to 'direction' over which particles are initially sampled 

#### Observations
The energy density is in general a function of position (3 variables), wavevector (2 unknowns, because norm is related to the fixed frequency) and time. For simplicity, the wavevector is assumed to only depend on the angle between the propagation direction and the position vector. The choice in this code has been to model the energy density as a 4-dimensional (3 dimensions for positions and one angle for direction) distribution evaluated at discrete times. To accelerate evaluation, the user can only consider two dimensions among the four ('x', 'y', 'z', 'directions'), while the other two are automatically integrated upon. The output energy density is therefore a 3-dimensional matrix (last dimension represents time).
The structure is composed of:
- 'x', 'y', 'z', 'directions' (Nx1 vectors): bins for space and propagation direction. In spherical coordinates (geometry.frame='spherical'), 'x', 'y', 'z' respresent respectively radius, azimuth (between -pi and pi) and elevation (between -pi/2 and pi/2). At least two of these must be 2x1 vectors to indicate over which values integration should be performed. The other two vectors should be Nx1 vectors (N>=2) indicating the boundaries of the bins over which distributions should be monitored.
- `time` (Mx1 vector): times at which distributions should be monitored.

#### Material
Material can be specified using routine MaterialClass.
