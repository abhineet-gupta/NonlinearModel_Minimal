1. calc_Tsg()
   This fuction provides the grid transformation matrix between the structural grid and the spline grid.

2. calc_Tks()
   It provides the grid transformation matrix between the spline grid and the aerodynamic grid based on a radial basis function for an infinite flat plate.


3. modal_aero()

  This function generates the generalized aerodynamic matrices for specified frequencies. 


4. rational_approx()
   This function fits the 8 GAM matrices in a least squares manner to a transfer function of the form:

  A(s) = A0 + s*A1 + (s^2)*A2 + state_space(A_lag,B_lag,C_lag) 

5. getDiff()
   This function computes the differentiation matrices which relate the positions and velocities of all the panels to the downwash experienced by them at their control points. 