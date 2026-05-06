# 🚀 Triaxial Accelerometer Static Calibration

This repository contains the MATLAB implementation and the comprehensive report for the **Optimal Estimation** course midterm assignment (Beihang University). 

The project focuses on the **self-contained static calibration** of a triaxial MEMS accelerometer, estimating deterministic errors (scale factors and zero-g biases) without relying on built-in MATLAB non-linear optimization functions.

## Overview
In real-world aerospace and industrial applications, sensors are affected by stochastic noise and deterministic errors. This project explores, implements, and rigorously compares two different calibration estimators:

1. **Linear Calibration (Attitude-Dependent):** Uses the Ordinary Least Squares (OLS) algorithm. It relies on a direct measurement model and requires precise knowledge of the sensor's attitude angles.
2. **Non-Linear Calibration (Attitude-Independent):** Uses a custom implementation of the **Gauss-Newton** algorithm. It relies on the inverse measurement model and exploits a physical invariant (the scalar norm of the local gravity vector $||A|| = g$) to decouple intrinsic sensor errors from spatial orientation.

## Key Features & Analysis
- **Custom Optimization Engine:** The Gauss-Newton algorithm and its Jacobian matrix were analytically derived and implemented from scratch.
- **Monte Carlo Sensitivity Analysis:** Extensive stress tests were conducted to evaluate algorithm robustness against:
  - Increasing stochastic noise (AWGN) and the *noise rectification effect*.
  - Scarcity or poor spatial distribution of attitude measurements (matrix ill-conditioning).
  - Extreme random perturbations in the Initial Guess (up to $\pm 100\%$ error).
- **Breaking Point Identification:** The analysis successfully identifies the mathematical limits of the unregularized Gauss-Newton approach (e.g., singular Jacobian when initial scale factors approach zero).
