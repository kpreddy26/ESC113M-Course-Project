# Computation of Liquid-Phase Composition and Temperature of a Vapour-Liquid Mixture

##  Overview

This project is developed as part of the ESC113: *Computer Methods for Engineers* course at IIT. The objective is to compute the liquid-phase composition and temperature of a binary vapour-liquid mixture at equilibrium, using classical numerical techniques and an interactive MATLAB-based GUI.

The application demonstrates the use of root-finding methods — Newton-Raphson, Bisection, and Regula Falsi — to solve thermodynamic equilibrium equations involving Raoult’s Law and Antoine’s Equation.

---

##  Problem Statement

Given:
- Two chemicals and their partial pressures at equilibrium
- Antoine coefficients for saturation pressure estimation

Find:
- Liquid-phase mole fractions (`x1`, `x2`)
- Temperature of the mixture at equilibrium

Assumptions:
- Ideal gas and ideal liquid behaviour
- Raoult’s Law holds

---

##  Numerical Methods Used

- **Newton-Raphson (Multivariable)**: For computing mole fraction (`x1`) and temperature simultaneously.
- **Newton’s Method (Single Variable)**: For refining temperature using a precomputed `x1`.
- **Bisection Method** and **Regula Falsi**: For comparison of root-finding performance.

---

##  App Features

- **User Inputs:**
  - Choice of two chemicals (from Acetone, Benzene, Toluene, Ethylbenzene)
  - Partial pressures for each chemical
  - Tolerance level
  - Initial temperature guess

- **Interactive Elements:**
  - Dropdown menus for chemical selection
  - Spinner for temperature
  - Radio buttons to select method
  - PLOT / CLEAR buttons for plotting convergence

- **Output:**
  - Estimated mole fractions using Newton-Raphson
  - Temperature computed via all three methods
  - Error convergence plots
  - Dew-point and bubble-point curves

---

##  Visualizations

- **Error Convergence Plots**: Iteration vs error for each method
- **Phase Curves**: Bubble-point and dew-point curve plotting at the calculated pressure and temperature

---

##  Technologies Used

- MATLAB 
- Numerical computation methods
- Antoine Equation
- Raoult’s Law

---

##  Reference

- Dorfman, K.D., & Daoutidis, P. (2017). *Numerical Methods with Chemical Engineering Applications*. Cambridge University Press.
- Finlayson, B. A. (2006). *Introduction to Chemical Engineering Computations*. John Wiley & Sons.

---

##  Future Improvements

- Support for non-ideal mixtures via activity coefficient models
- Import Antoine constants from databases like NIST
- Extension to multicomponent systems
- More chemicals and user-defined Antoine constants

---

## License

This project is developed for educational purposes only as part of the ESC113 course.

