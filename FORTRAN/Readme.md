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
 --------------------------------------------------------------------------------

This is a demonstration code for binary mixtures onlye, for more efficient
implementations please contact rafael.pelegrini [at] ufrgs.br.

 PLEASE CITE AS
   * Soares and Gerber (2013), Ind. Eng. Chem. Res. DOI:10.1021/ie400170a;
   * Soares et al. (2013), Ind. Eng. Chem. Res. DOI:10.1021/ie4013979.

The **src** folder contains FORTRAN source code.

The resulting executable can compute activity coefficients with the F-SAC model
for multi-component mixtures. For binary mixtures, in particular, a VLE diagram
can also be optionally computed.

The **pars** folder contains the model parameters for a limited number of
of substances.

More details for the included source files:

 - Main.f90 - the main program. This can be modified or replaced in case of using this model implementation in another program; 

 - Allocation.f90 - contains subroutines that allocate all the vectors and matrices, set the variables and compute the variables that do not change during program execution;

 - Activity.f90 - contains subroutines that compute the activity coefficient according to the F-SAC model for a given composition;

 - CompsData.f90 - reads groups and subgroups data contained in lib folder files;

 - RaoultPressure.f90 - computes the system's total pressure and vapour composition in equilibrium with the liquid phase (constant temperature). It uses modified Raoult's law. 

