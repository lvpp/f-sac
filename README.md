# F-SAC - Functional-Segment Activity Coefficient model

The F-SAC model is an activity coefficient model for mixtures. Its purpose is similar to the UNIFAC model.

The F-SAC model is based on the concept of functional groups but with the interaction energies between groups determined as in the COSMO-RS theory.

Most model parameters are for functional groups alone (not pairwise). This increases the predition power of the model and reduces the dependency from experimental data for model calibration. The only binary parameters of the model are the hydrogen bond (HB) formation energies.

This repository currently contains [FORTRAN](https://github.com/lvpp/f-sac/tree/master/FORTRAN) demonstrative code and a subset of parameters. For optimized code or updated parameter values please contact rafael.pelegrini [at] ufrgs.br. 

The relevant literature is:
 - [Functional-Segment Activity Coefficient model. 1. Model formulation](http://dx.doi.org/10.1021/ie400170a). Industrial & Engineering Chemistry Research, v. 52, p. 11159-11171, 2013.
 - [Functional-Segment Activity Coefficient model. 2. Associating mixtures](http://dx.doi.org/10.1021/ie4013979). Industrial & Engineering Chemistry Research, v. 52, p. 11172-11181, 2013.
 - [Simultaneous Correlation of Infinite Dilution Activity Coefficient,Vapor-Liquid, and Liquid-Liquid equilibrium data with F-SAC](http://dx.doi.org/10.1016/j.fluid.2013.11.040). Fluid Phase Equilibria, v. 364, p. 31-41, 2014.
 - [Mutual Solubilities of Hydrocarbon-Water Systems with F-SAC](http://dx.doi.org/10.1016/j.fluid.2014.10.026). Fluid Phase Equilibria, v. 25, p. 122-133, 2014.
 - [Including dispersive interactions in the F-SAC model](http://dx.doi.org/10.1016/j.fluid.2016.02.043). Fluid Phase Equilibria, v. 426, p. 56-64, 2016
 - [Prediction of water solubilities in hydrocarbons and oils using F-SAC coupled with SRK-EoS](http://dx.doi.org/10.1016/j.fluid.2016.08.001). Fluid Phase Equilibria, v. 427, p. 394-405, 2016.
 - [Extension of the F-SAC model to ionic liquids](https://doi.org/10.1016/j.fluid.2018.08.018). Fluid Phase Equilibria, v. 477, p. 87-97, 2018.
