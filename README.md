# Numerical-Simulation-of-LGRBs-Varying-Density-profile
This project performs numerical simulations of Long Gamma-Ray Burst (LGRB) relativistic jets propagating through varying density profiles of progenitor Wolf–Rayet stars. We utilize the **PLUTO code (v4.0)** to solve the equations of special relativistic hydrodynamics.

A key feature of this repository is the implementation of a **Black Hole (BH) spin evolution model**, based on the premise that relativistic jets extract angular momentum and "grind" the BH spin to a halt in Magnetically Arrested Disk (MAD) states.

## 🚀 Getting Started

### Prerequisites
* **Compiler:** `gcc` for serial runs or an MPI-wrapped compiler (e.g., `mpicc`) for parallel runs.
* **Libraries:** Standard C libraries. Use `ldd ./pluto` to check library links after compilation.
* **Python 3.x:** Required for the `setup.py` script and post-processing analysis.

### Project Structure
* `init.c`: Implements initial conditions, Wolf-Rayet density profiles, and boundary conditions.
* `definitions.h`: Preprocessor directives (physics module, geometry).
* `pluto.ini`: Runtime parameters (grid size, CFL, output frequency).
* `setup.py`: Python script used to configure the physics and makefile.
* `/Analysis`: Python scripts for data visualization and tracking BH spin evolution.
 