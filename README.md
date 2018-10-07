# Two-Energy-Gap-based-BTK-Theory
The matlab program is coded by github user Sora-Kasugano. All rights reserved.

This program is intended for the experimental fit of modified BTK theory considering the Two-Energy Gap assumption. Related work stemmed from the paper: PhysRevB.63.104510. The following codes are based from this paper.

Tips of this program:
1: ''V'' stands for a group of DC bias voltage data. By default, the unit of V is mV.
2: ''aa'' stands for a row vector, which contains the following fit parameters.
a= aa(1); b= aa(2); Delta1= aa(3); Delta2= aa(4); Gama=aa(5); Z= aa(6); P= aa(7); npanel= aa(8); T= aa(9);

     a ,b          = lower and upper limits of the integral.
     Delta         = Superconducting energy gap. In Two-Energy Gap based modified BTK theory, Delta1 is no more than Delta2 and                          their units are mV.
     Gama          = Inelastic scattering factor with its unit mV.
     Z             = Barrier factor
     P             = Spin Polarization
     T             = Experimental Temperature (K)
     npanel        = number of panels to use in the integration, with Total number of nodes = 2*panel + 1

The output form are [S,DI1], where S represents the DC bias voltage, and DI1 the normalized differential conductance.
