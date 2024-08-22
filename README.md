# Test Case Damping Identification
This repository contains test cases for various damping identification methods.

## Utilities

This folder contains helper functions.

## Linear SDOF Test Cases

This folder contains scripts for testing linear damping identification methods on time histories representative of a linear SDOF system without noise. 

The examined methods are the logarithmic decrement method, the matrix pencil method, and moving-block analysis. The methods are verified with the known exact solution for the system damping. 

The first two references below describe the logarithmic decrement and matrix pencil methods. The other reference describes moving-block analysis.

## References

Jacobson et al., "Evaluation of Time-domain Damping Identification Methods for Flutter-Constrained Optimization," Journal of Fluids and Structures, 2019. DOI: https://doi.org/10.1016/j.jfluidstructs.2019.03.011.

Kiviaho et al., "Flutter Boundary Identification from Time-Domain Simulations Using the Matrix Pencil Method," AIAA Journal, 2019. DOI: https://doi.org/10.2514/1.J058072.

Bousman and Winkler, "Application of the Moving Block Analysis," AIAA Paper 1981-653, 1981. DOI: https://doi.org/10.2514/6.1981-653.

## Contributors

Cristina Riso (Email: criso@gatech.edu)
