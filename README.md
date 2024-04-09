# Radiative Transfer Monte Carlo

## Introduction

This project proposes an approximation of solutions to the radiative transfer equation using the Monte Carlo method. It provides solutions for acoustic waves and elastic waves, both in 2D and 3D, which have been developed and implemented in this repository.

The library is primarily developed at the Laboratoire de Mécanique et d'Acoustique (LMA) in Marseille, France and Institut de Recherche en Génie Civil et Mécanique (GeM), Nantes University Nantes France.

- **Contact:** [Régis Cottereau](mailto:cottereau@lma.cnrs-mrs.fr)
- **Contributors (in order of first commit):** [Shahram Khazaie](mailto:Shahram.Khazaie@univ-nantes.fr)

The code has been written in MATLAB based on the following article:

- Ryzhik, L., Papanicolaou, G., & Keller, J. B. (1996). Transport equations for elastic and other waves in random media. *Wave Motion*, 24(4), 327-370. [DOI](https://doi.org/10.1016/S0165-2125(96)00021-2)

## Installation

The installation can be done using the `AddRTMCLib.m` script. This function will add the folder containing the Radiative Transfer Monte Carlo repository to the MATLAB path, allowing access to the repository functions from any work folder. 

There are two ways to add the library:

1. Run `AddRTMCLib(baseFolder)`, where `baseFolder` is the folder path of the repository.
2. Run `AddRTMCLib` directly from within the repository folder.

## Literature Comparison

Our code is compared to the results described in the following papers:

1. J. C. J. Paasschens. Solution of the time-dependent Boltzmann equation, *Phys. Rev. E* 56(1), pp. 1135-1141 (1997). [DOI](https://doi.org/10.1103/PhysRevE.56.1135)
2. M. Hoshiba. Simulation of multiple-scattered coda wave excitation based on the energy conservation law. *Phys. Earth Planet. Int.* 67, pp. 123-136 (1991). [DOI](https://doi.org/10.1016/0031-9201(91)90066-Q)
3. H. Nakahara, K. Yoshimoto. Radiative transfer of elastic waves in two-dimensional isotropic scattering media: semi-analytical approach for isotropic source radiation. *Earth Planets Space* 63, pp. 459-468 (2011). [DOI](https://doi.org/10.5047/eps.2011.03.006)
4. H. Sato. Multiple isotropic scattering model including P-S conversions for the seismogram envelope formation. *Geophys. J. Int* 117, pp. 487-494 (1994). [DOI](https://doi.org/10.1111/j.1365-246X.1994.tb03946.x)

The function `mainValidation(type)` can be used to run the comparison between our code and the literature presented above. The input of the function is: 'all', '2dIsotropicAcoustic', '3dIsotropicAcoustic' or '2dIsotropicElastic'.

## Usage

Two main examples are provided in this repository: `mainAcoustics.m` and `mainElastics.m`. These examples are tailored for solving the acoustic and elastic wave equations, respectively.

Additionally, there are two functions available for different problem scenarios:

radiativeTransferUnbounded: This function is suitable for full-space problems or scenarios without boundaries. It computes the radiative transfer using the specified input parameters and returns the results in a structure named obs.

Usage: `obs = radiativeTransferUnbounded(geometry.dimension, source, material, observation);`

radiativeTransferAcoustics: This function is designed for problems with boundaries. Similar to radiativeTransferUnbounded, it calculates the radiative transfer based on the provided input parameters and returns the results in the obs structure.

Usage:  `obs = radiativeTransferAcoustics(source, material, observation, geometry);`

Here is the format for the input data required by these functions:

### Input Data

The input is given by five MATLAB structures: `source`, `material`, `observation`, `geometry`, `Plotting`.

#### Observations

The structure is composed of:
- `dr`: size of bins in space
- `time`: observation times
- `Ndir`: number of bins for directions

#### Source
- `numberParticles`: number of particles in the simulation
- `position`: position of the source
- `lambda`: standard deviation of a Gaussian source

#### Material
- Acoustic Case:
  - `acoustics`: true
  - `v`: celerity of the homogeneous media
  - `Frequency`: frequency in Hz
  - `correlationLength`: Correlation length of the media fluctuation 
  - `spectralType`: Correlation of the media fluctuation 
  - `coefficients_of_variation`: defines the coefficients of variation of kappa (bulk modulus) and rho (density)
  - `correlation_coefficients`: defines the correlation coefficient between kappa (bulk modulus) and rho (density)

- Elastic Case:
  - `acoustics`: false
  - `vp`: pressure celerity of the homogeneous media
  - `vs`: shear celerity of the homogeneous media
  - `Frequency`: frequency in Hz
  - `correlationLength`: Correlation length of the media fluctuation 
  - `spectralType`: Correlation of the media fluctuation 
  - `coefficients_of_variation`: defines the coefficients of variation of lambda, mu (Lamé coefficients) and rho (density)
  - `correlation_coefficients`: defines the correlation coefficient between (lambda,mu), (lambda,rho), and (mu,rho)

#### Geometry

- Acoustic Case:
  - `type`: type of geometry. 'type' is either 'fullspace', 'halfspace', 'slab', or 'box'. 'halfspace' is defined for z<0. 'slab' is bounded between two planes of constant z. In 2D, 'box' is bounded in x and z. The elastic material only support
  - `size`: [width depth height]
  - `dimension`: dimension. It can be 2 or 3

- Elastic Case:
  - It is always full space
  - `dimension`: dimension. It can be 2 or 3

#### Plotting 
- equipartition: true - it plots the equipartition of energy
- movieTotalEnergy: true - it generates a gif
- movieDirectionalEnergy: true - it generates a gif
- sensors: positions in a matrix form with size of sensors by dimension