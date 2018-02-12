# spherocylinder
Metropolis Monte Carlo code for simulation of rigid spherocylinders with adhesive end-groups.
The code performs a metropolis monte carlo simulation in the NPT ensembe of a system of rigid cylinders with spherical end-caps (spherocylinders). The spherocylinders have sticky spots with angular specificity (angle theta) on their ends.
This project is intended as the starting point for custom simulations of the same or simular systems. For example, monte carlo move types can be added or removed, and the details of the simulation process can be adjusted.

## reqirements
Gnu Scientific Library (gsl >= 2.1)

Gnu Compiler Collection (gcc >= 6.1)

# usage
## compiling
compile using makefile

## running
This project is intended as the starting point for a custom system simulation. A sample implementation can be found in spcrun.cpp. Given an initial condition as a VTF file, it simulates the system by slowly changing the interaction strength beta and the system pressure from an initial value to a final value. The details of the spherocylinder geometry are given in the initial condition file, and the interaction angle theta is submitted as a parameter.
