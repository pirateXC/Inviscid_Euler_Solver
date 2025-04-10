# 2D Inviscid Euler Solver using Finite Volume Method (Steger-Warming Flux Splitting)

## Overview

This project implements a 2D finite volume solver for the inviscid Euler equations using the **Steger-Warming flux vector splitting** method. The goal is to numerically solve various compressible flow problems and validate the solver across three canonical test cases involving oblique shock waves and expansion fans.

The solver targets the conservative form of the Euler equations, discretized on a structured grid. The flux vector splitting technique enables stable upwinding in supersonic and transonic regimes by decomposing the fluxes based on the sign of the eigenvalues.

## Objectives

- Implement a conservative finite volume method for 2D Euler equations.
- Use Steger-Warming flux vector splitting for accurate shock capturing.
- Simulate and analyze:
  - Supersonic compression corners (oblique shocks)
  - Supersonic expansion corners (Prandtl-Meyer fans)
  - Oblique shock reflection from a wall
- Validate numerical results against analytical oblique shock and expansion theory.

## Governing Equations

The 2D Euler equations in conservative form are:
\[
\frac{\partial \mathbf{U}}{\partial t} + \frac{\partial \mathbf{F}}{\partial x} + \frac{\partial \mathbf{G}}{\partial y} = 0
\]
Where:
- \(\mathbf{U} = [\rho, \rho u, \rho v, \rho E]^T\)
- \(\mathbf{F}, \mathbf{G}\): Convective flux vectors in x and y

## Numerical Method

- **Spatial Discretization**: Finite Volume Method (FVM) on a Cartesian mesh
- **Flux Evaluation**: Steger-Warming flux vector splitting
- **Temporal Scheme**: Explicit time stepping (Runge-Kutta or forward Euler)

## Initial and Boundary Conditions

- **Initial Condition**: Uniform supersonic freestream conditions:
  - Mach number: \( M_\infty = 2.0 \)
  - \( \rho_\infty = 1.0 \), \( u_\infty = M_\infty \), \( v_\infty = 0 \)
  - \( p_\infty = 1.0 \)

- **Boundary Conditions**:
  - Inlet: Supersonic inflow (Dirichlet)
  - Outlet: Supersonic outflow (Extrapolation)
  - Top wall: Slip wall (Neumann)
  - Bottom wall: Slip wall or geometry-dependent boundary (e.g., wedge)

## Test Cases

### Case 1: Oblique Shock from Wedge

- Flow over a compression ramp (wedge with deflection angle).
- Expected: Oblique shock attached at the wedge corner.
- Analytically validated using oblique shock relations.

### Case 2: Expansion Corner (Prandtl-Meyer Fan)

- Flow over an expansion ramp.
- Expected: Prandtl-Meyer expansion fan emanating from the corner.
- Analytically validated using expansion theory.

### Case 3: Oblique Shock Reflection

- A reflected oblique shock interacting with a wall.
- Shock impinges and reflects, producing a complex flow pattern.
- Demonstrates capability of solver to capture shock-shock and shock-wall interactions.

## References

- [1] **AIAA 2001-2609**, Radespiel & Kroll: *Accurate flux vector splitting for the Euler equations*.
- [2] Anderson, J. D. *Modern Compressible Flow: With Historical Perspective*, 3rd ed.
- [3] Hirsch, C. *Numerical Computation of Internal and External Flows*, Vol. 2.
- [4] Lecture slides and class materials from AE 6042 (Spring 2025).

---
