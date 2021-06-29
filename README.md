# SWAT_SWC code and executables
#### See [Pignotti et al. (2021)](https://onlinelibrary.wiley.com/share/author/EARIJMERNHWYY57XECAJ?target=10.1002/hyp.14034) <link> for details about the theory and development behind this code.

## This repository contains three directories
1. compiled (modified) SWAT executables (both Windows and Linux)
2. modified source code (version 664) used to compile executables
3. example input files


### 1. Executables
Two compiled executables are made available for Windows (SWAT_nFC.exe) and Linux (swat664_nfc) systems in /gpignotti/swat_swc/tree/master/executables

Please note: For  comparison of output soil water values, the values in the output.swr file have been modified to print total soil water content, rather than available water content as specified by the default SWAT output. That is, the output.swr file now records available soil water plus wilting point water, not soil water content without wilting point water content.

Also note: Changes made to the code provide the user the ability to print volumetric soil water content for specified HRUs and associated layers in file.cio. HRU output variables 80-89 will print volumetric soil water content for soil layers 1-10 respectively for user-specified HRUs. Volumetric soil output is useful when comparing to observed soil moisture information which is generally recoded volumetrically as opposed to the equivalent water depth used in SWAT.

### 2. Modified source code
The main focus of modifying the SWAT code was to incorporate alternative equations to calculate soil water percolation based off the Campbell and van Genuchten approximations of hydraulic conductivity as described in the body of the paper. Additionally, this necessitated: 1) calculating the parameters for the equations based on soil properties, 2) creating an hourly loop for percolation, 3) better constraining maximum and minimum percolation and soil water content, and 4) printing new output files.

Please see the document: [Code_Edits_Summary.pdf] for a full description of all SWAT source code files edited.


### 3. Example input files
There are two new flags that can be used to specify which soil water equation to use and to turn on printing of extra variables.

(1) In the basins.bsn file an two extra lines are needed, one on line 17 and one at the bottom line of the file to specify percentage change to the exponential value of the Campbell and van Genuchten equations and the soil water percolation equation used by SWAT as follows:

Line 17:

0.000    | SWCEXP : Multiplicative factor used for calibration of b/m parameter in CA and VG equations

SWCEXP: A flag value of zero will not change the exponent value and will simply use the default value calculated by the PTFs as specified from the IPERC flag. Non-zero values will apply a percentage change to the exponent, either positive or negative as specified by the user.

Last line:

1	| IPERC: Changes the percolation method (0=default, 1=CA-RA, 2=VG-RA, 3=CA-CO, 4=VG-CO) 

IPERC: A flag value of 0 will use the default SWAT soil percolation equation. This is not recommended, rather the user should simply run an officially released SWAT executable if the default is preferred. Flag values of 1, 2, 3, and 4 indicate the Campbell-Rawls, van Genuchten-Rawls, Campbell-Cosby, and van Genuchten-Cosby equations respectively. Although, other options have been coded (5-8), users are cautioned against their use as they have not undergone testing to date.

(2) The second flag is in file.cio and turns on and off printing of additional soil water variables and is similarly added as a line at the bottom of the file as follows:

1	| IWRITEPERC: Code for printing new output soil files

Where 0 will not print the variables and 1 will print files for layer-based lateral flow, percolation, depth, and infiltration precipitation.

