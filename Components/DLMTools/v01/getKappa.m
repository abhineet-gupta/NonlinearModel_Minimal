function kappa = getKappa(x0,y0,z0,cosGamma,sinGamma,k,M)
%% Function to calculate kappa
% this function calculates kappa as defined on page 3 of reference 1. The
% All the formulae are from page 1 of the reference.
% kappa = (r1^2) * K where K is the incremental Kernel function
% K = (K1T1 + K2T2)*exp(-jwx0/U)/(r1^2), where w is oscillation frequency
% variables passed to the function:

% x0 = x - chi (x is location of collocation pt (3/4th chord pt), chi is
%               location of doublet)
% y0 = y - eta
% z0 = z - zeta
% cosGamma: cosine of panel dihedral
% sinGamma: sine of panel dihedral
% k = w/U,  U is freestram velocity
% M: Mach no.

% Reference papers
% ref 1: Albano and Rodden - A Doublet-Lattic Method for Calculating 
%        Lift Distributions on Oscillating Surfaces in Subsonic Flows
%
% ref 2: Watkins, C. E., Hunyan, H. L., and Cunningham, H. J., "A Systematic 
%        Kernel Function Procedure for Determining Aerodynamic Forces on Oscillating 
%        or Steady Finite Wings at Subsonic Speeds," R-48, 1959, NASA.
%
% ref 3: Blair, Max. A Compilation of the mathematics leading to the doublet 
%        lattice method. No. WL-TR-92-3028. WRIGHT LAB WRIGHT-PATTERSON AFB OH, 1992.
%
% Functions called: getI1(u1,k1)
%                   getI2(u1,k1) 
% Code by:
% Aditya Kotikalpudi
% University of Minnesota
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% declare all variables as defined in reference 1, page 1
%z0 = zeros(size(y0));
r1 = sqrt((y0.^2) + (z0.^2));
beta2 = (1-(M^2));
R = sqrt((x0.^2) + (beta2*(r1.^2)));
u1 = ((M*R) - x0)./(beta2*r1);
k1 = k*r1;
j = sqrt(-1);

%% calculate T1, T2


cos1 = repmat(cosGamma,1,length(cosGamma));
cos2 = repmat(cosGamma',length(cosGamma),1);
sin1 = repmat(sinGamma,1,length(sinGamma));
sin2 = repmat(sinGamma',length(sinGamma),1);

T1 = cos1.*cos2 + sin1.*sin2;
T2 =  (z0.*cos1 - y0.*sin1).*(z0.*cos2 - y0.*sin2)./(r1.^2);

%% calculate I1 & I2
I1 = getI1(u1,k1);
I2 = getI2(u1,k1);

%% get kappa_temp
kappa_temp1 = I1 + ((M*r1).*exp(-j*(k1.*u1))./(R.*sqrt(1+(u1.^2))));
kappa_temp2 = -3*I2 - (j*k1.*(M.^2).*(r1.^2).*exp(-j*k1.*u1)./((R.^2).*sqrt(1+u1.^2))) - ...
              (M*r1.*((1+u1.^2).*((beta2*r1.^2)./R.^2) + 2 + (M*r1.*u1./R))).*exp(-j.*k1.*u1)...
              ./(((1+u1.^2).^(3/2)).*R);
          
kappa_temp = kappa_temp1.*T1 + kappa_temp2.*T2;

%% Resolve the singularity arising when r1 = 0, ref 2, page 7, Eq 18
rInd = find(r1==0);  
kappa_temp(rInd(logical(x0(rInd)<0))) = 0;
kappa_temp(rInd(logical(x0(rInd)>=0))) = 2;

kappa = kappa_temp.*exp(-j*k*x0);


return
