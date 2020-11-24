This folder contains scripts to demonstrate the use of the DLM software on the BFF aircraft. 

The folders are:
Data_files: It contains all the .mat files required to construct the aerodynamic model-
1. AeroGridData_v01.mat: contains panels information required to solve for DLM solution.

2. GAM_data_v01: contains transformation matrices which create the final generalized                  aerodynamic matrices in modal coordinates 

3. StrucDynData_v01: contains data from the FE model

4. AircraftData_v01: contains aircraft properties like mass, inertia and reference                                  geometries

The scripts in the folder are:

AICgen: The script which gives the final .mat file to be plugged into the simulation, uses all         the data contained in the Data_files folder and scripts in DLM_code and GAM_generation         folders to get the final outputs.