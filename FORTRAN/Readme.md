--------------------------------------------------------------------------------
 Copyright (c) 2011-2013, Federal University of Rio Grande do Sul
 All rights reserved.

 This software is subject to a BSD License, please see the License.txt
 file for more information.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED.

 Authors: Luiz Felipe Kusler Possani
          Rafael de Pelegrini Soares

 This is a demonstration code, for more efficient implementations
 please contact rafael@enq.ufrgs.br.

 PLEASE CITE AS
   * Soares and Gerber (2013), Ind. Eng. Chem. Res. DOI:10.1021/ie400170a;
   * Soares et al. (2013), Ind. Eng. Chem. Res. DOI:10.1021/ie4013979.
--------------------------------------------------------------------------------

This folder contains:

 - License.txt - license file;

 - FSAC1.exe - program executable. Insert the number of components in mixture, their respective names, system's temperature (constant) and composition. The program computes components' activity coefficients (gamma) in mixture. If a binary mixture is inserted, the program computes gamma profile by varying one component molar fraction from 0 to 1. Also, if component's vapour pressures are inserted, it computes the system's total pressure and vapour composition in equilibrium with the liquid phase;

 - lib - contains the following files:

  * FSAC1-comps.csv - contains name, the subgroups which the molecule is divided and other information; 
 
  * FSAC1-subgroups.csv - contains the group that the subgroup belongs and the subgroup's volume and superficial area;

  * FSAC1-groups.csv - contains the positive and negative areas and information about what fraction of superficial area makes hydrogen bonds;
  
  * FSAC1-wHB.csv - contains energy of hydrogen bonding for pairs of groups.

The following files are created during program execution (only for binary mixtures):

 - Gamma.dat - gamma profile;

 - PressureProfile.csv - pressure profile (if the vapour pressures are inserted - if not, this file is not created).

Source files:

 - Main.f90 - the main program. This can be modified or replaced in case of using this model implementation in another program; 

 - Allocation.f90 - contains subroutines that allocate all the vectors and matrices, set the variables and compute the variables that do not change during program execution;

 - Activity.f90 - contains subroutines that compute the activity coefficient according to the F-SAC model for a given composition;

 - CompsData.f90 - reads groups and subgroups data contained in lib folder files;

 - RaoultPressure.f90 - computes the system's total pressure and vapour composition in equilibrium with the liquid phase (constant temperature). It uses modified Raoult's law. 

