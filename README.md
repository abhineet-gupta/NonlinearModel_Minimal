#Nonlinear Model

The repository contains a nonlinear simulation model of the mAEWing1 series of aircraft developed at University of Minnesota as part of the 'Performance Adaptive Aeroelastic Wing' project.
Linear models can be obtained by linearizing the model at speeds by running 'Simulation_6dof/NL_Sim/save_models.m'.
The scripts save the linearized models in the the structure named 'simulink_6dofVLM'.

Refer to the following thesis for further details on the simulation model.
```
@phdthesis{gupta2019flight,
  title={Flight Dynamics Model of a Small Flexible Aircraft},
  author={Gupta, Abhineet},
  year={2019},
  school={University of Minnesota}
}
```
