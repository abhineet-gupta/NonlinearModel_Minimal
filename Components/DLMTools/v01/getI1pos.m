function I1pos = getI1pos(u1,k1)
%% Function to get I1 for positive u1 values
% ref 1: Albano and Rodden - A Doublet-Lattic Method for Calculating 
%        Lift Distributions on Oscillating Surfaces in Subsonic Flows

% I1 described in eqn 3 of page1 of reference 1. Approximated as shown in
% page 3 of the reference. This implementation is only valid for u1>0. For
% u1<0, this function is still used to obtain I1 in an indirect manner as
% described in getI1.m, in which this function is called.
%
% Code by
% Aditya Kotikalpudi
% University of Minnesota
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% constants used in approximation

a1 = 0.101;
a2 = 0.899;
a3 = 0.09480933;
b1 = 0.329;
b2 = 1.4067;
b3 = 2.90;

%% solve for I1
j = sqrt(-1);
i1 = (a1*exp((-b1-(j*k1)).*u1)./(b1+(j*k1))) + (a2*exp((-b2-(j*k1)).*u1)./(b2+(j*k1)));
i2 = (a3./(((b3 + (j*k1)).^2) + (pi^2))).*(((b3+(j*k1)).*sin(pi.*u1)) + (pi*cos(pi.*u1))).*exp((-b3-(j*k1)).*u1);
I1_temp = i1 + i2;
I1pos = ((1-(u1./sqrt(1+u1.^2))).*exp(-j*k1.*u1)) - (j*k1.*I1_temp);