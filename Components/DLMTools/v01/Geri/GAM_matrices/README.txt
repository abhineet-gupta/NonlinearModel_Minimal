This folder contains scripts required in the GAM generation for the BFF aircraft.

1. spline_grid: It generates a spline grid based on the structural grid of the aircraft, which                 serves as the intermediate grid between the structural and aerodynamic model                 grid points. 

2. GAM_matrices: It generates the GAM_data.mat file using the scripts in the folder                             'GAM_generation' as well as the spline grid created by spline_grid.m