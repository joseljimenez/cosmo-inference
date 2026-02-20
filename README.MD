# Cosmo-Inference

Cosmo-Inference is a lightweight Bayesian framework for cosmological parameter estimation.  
It constrains background cosmological models using observational data through MCMC sampling.

Currently implemented models:
- ΛCDM
- Generalized Chaplygin Gas (GCG)

---

## Scientific Overview

The code performs parameter inference by fitting theoretical predictions of the Hubble expansion rate to observational datasets.

### ΛCDM Model

The Friedmann equation is given by:

H²(z) = H0² [ Ωm (1+z)³ + ΩΛ ]

Free parameters:
- H0
- Ωm

---

### Generalized Chaplygin Gas (GCG)

Equation of state:

p = -A / ρ^α

Energy density evolution:

ρ_GCG(z) = ρ0 [ As + (1 - As)(1+z)^{3(1+α)} ]^{1/(1+α)}

Free parameters:
- H0
- As
- α

ΛCDM is recovered in the limit α → 0.

---

## Observational Data

The framework supports background cosmology datasets such as:

- Type Ia Supernovae
- Cosmic Chronometers H(z)
- BAO measurements

Each dataset contributes through a Gaussian likelihood.

---

## Statistical Method

Parameter estimation is performed via Bayesian inference using MCMC sampling.

- Gaussian likelihoods
- Flat priors on cosmological parameters
- Modular sampler implementation

---
