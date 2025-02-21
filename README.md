# 1D Eulerian phytoplankton diversity model Readme

This file explains how to run the Eulerian phytoplankton diversity model expanded from a Nutrient-Phytoplankton-Zooplankton-Detritus (NPZD) model in a one-dimensional (1D) vertical water column. 
This model aims to be a Eulerian version of the individual-based Lagrangian module of the **IBM** repository.  

The model output has been compared and validated against the observations at the Bermuda Atlantic Time-series Study (BATS) station.

## Authors
Bingzhang Chen

## Contact email
[bingzhang.chen\@strath.ac.uk](mailto:bingzhang.chen@strath.ac.uk){.email}

## Brief summary of the study
This 1D Eulerian NPZD model was developed to compare the results with the 1D hybrid Eulerian-Lagrangian model developed in the **IBM** repository. 

Right now two modules have been developed: a simple NPZD model and a model with multiple size classes of phytoplankton and zooplankton. Both models employ the variable stoichiometry model of Geider et al. (1998) in which phytoplankton N:C and Chl:C ratios vary with the environment.

The size structure model largely follows Ward et al. (2012), but neglects the variable N:C ratio in zooplankton.

Currently, the model is configured to simulate the Bermuda Atlantic Time-series Study (BATS) station (64.17 $^\circ$ W - 31.67 $^\circ$  N) in the subtropical north-western Atlantic Ocean.

## Contributor

Bingzhang Chen is responsible for writing the Fortran code.

## LICENSE

All the data and codes are covered by the MIT license. Please see the LICENSE file for details.

# Metadata

## Software and packages

The codes are written in Fortran 90. The model codes have been tested in a x86_64 Red Hat linux system using intel fortran or gfortran. The code requires a number of input files (i.e., forcing) to run (see below).

## How to run the code

1. First, it should be checked thar ifort or gfortran compilers are installed on the working machine. To verify the installation, the user can type "ifort -v" or "gfortran -v" in the Terminal window.
	The installation instructions for ifort can be found on https://software.intel.com/content/www/us/en/develop/documentation/get-started-with-fortran-compiler/top.html.
	Those for gfortran can be found on https://fortran-lang.org/learn/os_setup/install_gfortran.

2. On the Terminal window, go to the directory where the model will run (we assume that the root directory is under home directory: ~/).

3. To download the code type: "git clone https://github.com/BingzhangChen/Euler1D.git".

4. To run the simple NPZD model, go to the working directory by typing "cd Euler1D/Run_simple".

5. To change the settings for the model run the file "job-Archie" must be configured. The file "job-Archie" is specifically used on the HPC cluster ARCHIE-WeSt of University of Strathclyde. In this file the user can define: (i) "Test = 0" for a fast run, usually for a formal model run for a large number of iterations; or (ii) "Test = 1" to run the model for debugging mode, which is much slower. However, for a user who does not use ARCHIE-WeSt, the user should select the right Fortran compiler to use, and also modify the compiler flags depending on the purpose in the script. The user also needs to correctly specify the NETCDF and openmpi directory. 

6. Edit the fortran file "params.f90" to specify the parameter values.

7. To compile the model type "./job-Archie", and an executable (**Euler1D**) will be generated.

7. Before running the model, there are two namelist files that the user needs to check and configure. In the file "time.nml" the user can define the time paramters for the model run (more details below). In the file "param.nml" are defined several plankton model parameters (more details below).

8. After defining all the model settings type the preferred command to run the model:

		“./Euler1D” to simply run the model with only one CPU; some details of the running will be shown on the screen.

		“./Euler1D \> out” to run the model with only one CPU and save the running details in the “out” file.

		“./Euler1D \> out & disown” to run the model with only one CPU in the background and save the running details in the “out” file.
	
	The model will generate a netcdf file with the Eulerian fields named "Euler.nc".

## Source codes

The following files are located in the directory Euler1D/src/:

- **Advection_center.f90**: subroutine to run advection of all Eulerian fields in the model. Here are also defined the boundary conditions.

- **Calc_PAR.f90**: subroutine to calculate vertical light attenuation based on chlorophyll profiles and attenuation coefficients.

- **Diff_center.f90**: subroutine to run diffusion of all Eulerian fields in the model.

- **forcing.f90**: module file containing several subroutines providing the external forcing (temperature, light, vertical eddy diffusivity) as a function of time.

	In this module, the subroutine *VERTICAL_LIGHT* computes the Photosynthetically Active Radiation (PAR) below the ocean surface (Anderson et al., 2015) based on the station latitude and the day of the year. Once calculated the PAR at the surface, it calls the *Calc_PAR* subroutine to compute the vertical light attenuation.

	The subroutine *extract_Kv* extracts profiles of vertical eddy diffusivity from external files (data from Vallina et al. (2017)).

	The subroutine *extract_WOAtemp* extracts temperature data from external forcing files (data source: WOA13).

- **Geider_Euler.f90**: this file contains all the main subroutines necessary to model the biological components of the Eulerian NPZD model. The subroutine *BIOLOGY* computes de Eulerian fields for nitrogen, zooplankton and detritus. This subroutine calls the model approach (Model_ID) selected by the user.

- **GMK98_simple.f90** (Model_ID = 1) computes the changes of phytoplankton cellular carbon, nitrogen and chlorophyll content following Geider et al. (1998).

- **GMK98_Size.f90** (Model_ID = 2) computes the changes of phytoplankton cellular carbon, nitrogen and chlorophyll content (Geider et al., 1998), considering the allometric relationships between the cell volume and the minimum cellular N:C ratio following Ward et al. (2012), and the allometric relationship between the cellular volume and the half-saturation constant for nutrient uptake following Edwards et al. (2012).
	
	
	
	
	- The subroutine *Par2PHY* calculates the total concentrations of phytoplankton carbon, nitrogen, and chlorophyll at each grid depth based on the super-individuals present.

- **grid.f90**: module file defining the model grid, which contains the *setup_grid* subroutine. Here the user can change the maximal depth and the vertical resolution.

- **gridinterp.f90**: a utility subroutine in which observational data, which might be given on an arbitrary, but structured grid, are linearly interpolated and extrapolated to the actual model grid.

- **initialization.f90**: subroutine to initialize all model variables (time settings, model parameters, nitrogen, zooplankton, detritus sinking rate, number of super-individuals, and nitrogen, carbon and chlorophyll cellular content of the phytoplankton cells), and forcing environments (temperature and diffusivity). This subroutine also creates the output nc files.

- **netcdf_IO.f90**: module file containing several subroutines to save model outputs to external netcdf files (Eulerian and Lagrangian files).

- **Main.f90**: main program running the model and timing the model run.

- **Makefile**: for compiling the Fortran codes and generating the executable **Euler1D**. If the user needs to generate new source files, they have to be added here.

- **../Run/params.f90**: module file declaring and assigning plankton model parameters. It initializes the parameters defined on the Run/param.nml file, and also defines the activation energy for phytoplankton growth and for zooplankton grazing, the maximal Chl:N ratio for phytoplankton, the value of rhochl before last sunset, and the mutation rate parameters mentioned above.

- **Readcsv.f90**: subroutine that reads the external forcing data files to compute the temporal variability of temperature and vertical diffusivity. These files ($*$.dat) can be located at the IBM/Run directory.

- **time_interp.f90**: subroutine that interpolates the external forcing fields from time series observations to model time.

- **Trait_functions.f90**: module file that contains several functions that calculate phytoplankton physiological rates from environmental factors (nitrogen, temperature and light) and traits (size, optimal temperature, and light). Here are also defined functions to estimate cellular carbon content from cellular volume, and vice versa, zooplankton prey palatability, and the photoinhibition.

- **Time_settings.f90**: module file declaring time-related variables. This module contains the subroutine *update_time* that converts each time step to the current second, second of the day, current DOY, etc.

- **timestep.f90**: this file contains the subroutines *TIMESTEP* and *UPDATE_PARTICLE_FORCING*, the main subroutines that update the status of the Eulerian fields and the Lagrangian particles (super-individuals), as well as the advection and diffusion at each time step.

- **tridiagonal.f90**: subroutine used in solving diffusion equations.

- **variables.f90**: module file declaring the state variables of Eulerian fields and the Lagrangian particles (super-individuals). In this file are defined the minimum and maximum zooplankton sizes, the total number of super-individuals, the initial number of phytoplankton cells per super-individual, and the mutation parameters.

	This module also contains the subroutine *UPDATE_PHYTO*, that updates the phytoplankton nitrogen, carbon and chlorophyll concentration by depth.


## Input data

The following files are located in the directory IBM/Run/:

- **BATS_Kv_time.dat**: this file contains the time file for the vertical eddy diffusivity forcing data obtained at the Bermuda Atlantic Time-series Study (BATS) station from Vallina et al. (2017).

- **BATS_Kv.dat**: this file contains the profiles of vertical eddy diffusivity data obtained at BATS station at each time step corresponding to BATS_Kv_time.dat.

- **BATS_temp_time.dat**: this file contains the time file for the temperature forcing data obtained at BATS station, extracted from World Ocean Atlas (WOA) 2013. 

- **BATS_temp.dat**: this file contains the temperature forcing data obtained at BATS station at each time point corresponding to BATS_temp_time.dat, extracted from WOA13.

- **BATS_NO3_Jan.dat**: this file contains the vertical profile of nitrate plus nitrite concentration at BATS station in January, extracted from WOA13. This file is used for initializing nutrient data of the model.

## Observational data

- **BATS_Primary_Production.csv**: this file contains the measurements of primary production at BATS station.

- **bats_NO3.csv**: this file contains the measurements of nitrate at BATS station.

- **bats_pigments.csv**: this file contains the measurements of Chlorophyll a at BATS station.

## Bash file to compile the code

- **job-Archie**: bash script to compile the Fortran files and generate the executable (./IBM). This script only works on Archie-West, the supercomputer cluster hosted by the University of Strathclyde. The user needs to modify the script to specify the path to mpi, the fortran compiler, and the netcdf library of their own computing facility.

## Namelist files for defining model parameters

- **param.nml**: namelist file containing several model parameters: maximal growth rate, the initial slope of the photosynthesis-irradiance curve, the half-saturation constant of nitrogen uptake, maximal zooplankton grazing rate, half-saturation constant of zooplankton grazing, zooplankton linear mortality, zooplankton gross growth efficiency, the fraction of unassimilated food ingested by zooplankton, the standard deviation of zooplankton feeding preference, the conversion rate from detritus to dissolved inorganic nitrogen and the detritus sinking. Here, the user can also define the Model_ID to select del modelling approach of interest.

- **time.nml**: namelist file controlling the parameters for model run. Here, the user can define the time parameters that establish the total number of simulation days (*NDay_Run*, in days), the time step (*dtsec*, in seconds), and the frequency at which the model outputs will be saved (*nsave*). For example, if the user wants the model outputs to be saved at a daily interval (i.e., every 86400 s), *nsave* should be equal to 86400/*dtsec*.

## R scripts to generate figures in the manuscript


- **FIG02_PhysicalForcing.m**: Matlab script to plot physical forcing (vertical eddy diffusivity, temperature, and light) driving the 1D model (Fig. 2 in the manuscript).

## Funding

This work is funded by a Leverhulme Trust Research, UK Project Grant (RPG-2020-389).

## References

Anderson, T. R., Gentleman, W. C., Yool, A. (2015) EMPOWER-1.0: An Efficient Model of Planktonic ecOsys-tems WrittEn in R. Geosci. Model Dev.8: 2231--2262. <doi:10.5194/gmd-8-2231-2015A>

Edwards, K.F., Thomas, M.K., Klausmeier, C.A., Litchman, E. (2012) Allometric scaling and taxonomic variation in nutrient utilization traits and maximum growth rate of phytoplankton. Limnology and Oceanography 57, 554--566. <doi:10.4319/lo.2012.57.2.0554>.

Geider, R.J., Maclntyre, H.L., Kana, T.M. (1998) A dynamic regulatory model of phytoplanktonic acclimation to light, nutrients, and temperature. Limnology and Oceanography 43, 679--694. <doi:10.4319/lo.1998.43.4.0679>.

Vallina, S. M., P. Cermeno, S. Dutkiewicz, M. Loreau, and J. M. Montoya. (2017) Phytoplankton functional diversity increases ecosystem productivity and stability. Ecol. Mod. 361:184–196.