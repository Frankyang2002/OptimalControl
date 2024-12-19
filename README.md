# Dodo Alive: Simulation and Control
This repository is part of the "Dodo Alive! - Resurrecting the Dodo with Robotics and AI: Simulation & Control" practical course

## Repository Structure
The following folders and files exist:
- **control**:   contains MATLAB scripts used for the impedance (joint/cartesian/polar) control and the hopping control

- **dynamics**:   contains the MATLAB script to derive the monopod dynamics ("*derive_dynamics.m*") as well as the derived dynamics functions and other useful functions, e.g. foot position, Jacobians, etc.

- **guards**:   contains the guard functions to detect events for discrete state transitions of the dynamics

- **img**:   contains generated images for the final report
    
- **parameters**:  contains the MATLAB script to set the monopod parameters, e.g. leg length, mass, etc. as well as the monopod parameters ("*monopod_parameters.mat*")

- **plots**:   contains MATLAB scripts to visualize the monopod dynamics simulation
    
- **resources**:   necessary for the simulink project

- **simulink**:   contains the simulink models of the monopod dynamics ("*monopod_model.slx*") and the impedance and hopping controllers ("*monopod_controller.slx*"), as well as the data dictionaries for both models to set parameters. The simulink model can later be used to export the MATLAB code to an embedded device

- **Monopod.prj**: MATLAB project used to organize the code

- **simulate_monopod.m**: main script to simulate the monopod dynamics and controllers

## How To Use the Code
To use the code for the simulation of the monopod hopper, first start the simulink project "*Monopod.prj*" using MATLAB.
The project automatically adds the needed paths for the code to work.
After opening the project, a simulation can be started by running the script "*simulate_monopod.m*".
To change the monopod parameters use the script "*set_monopod_parameters.m*" in the parameters folder. Afterwards, you have to run the script "*derive_dynamics.m*" inside the dynamics folder to generate new monopod dynamics functions.
