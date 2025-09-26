# Pinns for Economics

The aim of this repository is to translate few standard PDE problems of economics and finance (HJB; Kolmogorov-forward/Fokker-Planck) to Physics-informed neural networks.
There will be a total of 6 Matlab examples.

## Content of Repository

* [Introductory lecture on Pinns](Pinn_lecture)
* [6 benchmark models, implemented in Matlab](matlab_examples)
    - [01: Partial equilibrium model](matlab_examples/01_partial_equilibrium_diffusion), explained in this [note here](matlab_examples/notes.pdf)
    - [02: Partial equilibrium model with discrete choices](matlab_examples/02_partial_equilibrium_diffusion), explained in this [note here](matlab_examples/notes.pdf)
    - [03: General equilibrium model; HJB plus Kolmogorov Forward](matlab_examples/03_general_equilibrium), explained in this [note here](https://benjaminmoll.com/wp-content/uploads/2020/02/HACT_Numerical_Appendix.pdf). Honestly Ben’s numerical appendix I sent above is perfect documentation of model 3. Read sections 1-3


## Aims

* Replicate all the results from the Matlab benchmark with high accuracy with PINNs.
* Carfully check what boundary conditions work fine (hard/soft boundary conditions).
* Also look at [BSDE solver methods](https://arxiv.org/abs/2505.17032), and solve models also with this deep learing approach.
* Look at this website: 
* Keep in mind that we want to build up know-how so we can add dimensions later to those models; boundary conditinons should not suffer from the curse of dimensinality.
* **IMPORTANT: we are interested accurate policies, that is, the slope of the function of the problems.** 

    
