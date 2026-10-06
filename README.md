# Phase-field-frictional-heating-vumat-vumatht

Coupled VUMAT–VUMATHT subroutines for phase field fracture with frictional and crack-propagation
heat sources in ABAQUS/Explicit, supporting Miehe and hybrid strain-energy splits.

## Overview

This repository provides an ABAQUS/Explicit implementation of the phase field fracture method,
extended with a frictional heating source along sliding crack faces. The phase field is
represented via the temperature degree of freedom in a Coupled temperature-displacement analysis
(the thermal analogy), solved through two independent user subroutines:

- **VUMAT** — solves the mechanical (displacement) problem: strain-energy decomposition, degraded
  stress update, crack driving force $\mathcal H$, and frictional heat generation rate.
- **VUMATHT** — solves the phase-field (temperature-analogy) problem: evolves $\phi$ via the
  heat-conduction analogy and feeds $\nabla\phi$ back to the VUMAT for the next increment.

The two subroutines communicate through a shared Fortran module (`pf_shared`), staggered by one
increment, consistent with the explicit time integration scheme.

Two variants are provided, differing only in the strain-energy split used to drive the crack
driving force:

| File | Strain-energy split |
|---|---|
| `VUMAT_VUMATHT_Miehe.f` | Miehe-type spectral split (reference/comparison model) |
| `VUMAT_VUMATHT_Hybrid.f` | Volumetric + principal deviatoric (hybrid) split |

Full theory, implementation details, SDV layout, and usage instructions are in
[`documentation.pdf`](documentation.pdf).

## Repository structure
