The code has various components related to unsteady aerodynamics, steady aerodynamics, gridding the geometry and code verification. All the components except for gridding can be used individually and independently.

The various folders are:

2. GAM_generation: It contains the functions used to generate transformation matrices which                    convert DLM solutions at discrete frequencies into a transfer function                    form to be plugged into the nonlinear simulation.

3. verification: contains a basic validation test of the DLM software

4. BFF_example: contains an example based on the BFF aircraft to demonstrate the use of the                 software, generates the .mat files required for the simulation

 
Description of the function of the core DLM implementation:

1. getAIC()
  The getAIC() function computes the AIC matrix for a given reduced frequency. It computes the downwash matrices from the DLM (getDLM()) for unsteady effects, after subtracting out the steady part of the solution (zero reduced frequency). The steady aerodynamics solution is added separately via a VLM implementation (getVLM()). 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

2. getVLM()

  The VLM function can be used independently for any aircraft, provided the inputs are in the compatible form. A brief description of inputs to the VLM code are provided:
 
  getVLM() inputs:

a. P0, P1, P3: These are the coordinates of the collocation/control point, southern end of vortex line on quarterchord, and northern end respectively.

b. PAreas: vector containing panel areas.

c. M: Mach no. (not used here, but can be)

d. n_hat_w: Vector of cosine of panel dihedral angle. It's value is 1 for all wing panels (zero dihedral angle) and 0 for winglet panels. 

e. n_hat_wl: Vector of sine of panel dihedral angles. It's 1 for all winglet panels (which are at 90 degrees to the wing) and zero for wing panels. 

  The output from the function is the downwash matrix D, which can be inverted to obtain the steady aerodynamics AIC matrix. It has to be remembered that the VLM enforces a zero normal flow boundary condition on each panel and computes the induced downwashes perpendicular to panel surfaces. This means that when the resulting AIC matrix is used for lift computation, it has to be used carefully, as follows-

                           p = [AIC]*w   (1)

where 'p' is row vector of pressure difference across all panels, 'w' is upwash angle experienced by each panel. It has to be remembered that 'w' represents the normal flow experienced by each panel of the aircraft. Therefore, it depends on the angle of attack (alpha), sideslip angle (beta) and the panel dihedral angle. The expression for 'w' will be:

                 w = (alpha*n_hat_w) + (beta*n_hat_wl)   (2)

  So, for example, if the BFF aircraft is at an angle of attack alpha (beta = 0), only the wing panels will experience an upwash. When the BFF aircraft is under sideslip beta, with alpha = 0, only the winglets experience the upwash. 

  After computing the RHS of equation 1, the resulting pressure difference for i'th panel can be used to obtain lift across that panel using equation 3
                      
                     L(i) = p(i)*qbar*Area(i)    (3)

  Again, it has to be remembered that the row vector L obtained here is prependicular force acting on each panel. To obtain lift at each panel along the global z axis,
                        
lift force:             Lift = L*n_hat_w

side force:                Y = L*n_hat*wl      (4)

The resulting vectors can then be summed up to obtain total forces along z & y axes.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

3. getDLM()
   The getDLM() function is used to obtain the Downwash matrix for unsteady aerodynamics, after the steady solution has been subtracted from it. The resulting downwash matrix is a function of the reduced frequency, typically defined as 
                            
                           k = cw/2V    (5)
where 'w' is the oscillating frequency of the wing, 'c' is the reference chord and 'V' is the airspeed. However, in the code, the reduced frequency is specified as
                           
                           k = w/V     (6)

Therefore, this factor has to be kept in mind while specifying a reduced frequency value for an aircraft in conventional sense as in equation 5.
