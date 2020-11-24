function I2 = getI2(u1,k1)
%% Function to get I2
% ref 1: Albano and Rodden - A Doublet-Lattic Method for Calculating 
%        Lift Distributions on Oscillating Surfaces in Subsonic Flows
%
%ref 3: Blair, Max. A Compilation of the mathematics leading to the doublet 
%       lattice method. No. WL-TR-92-3028. WRIGHT LAB WRIGHT-PATTERSON AFB OH, 1992.
%
% I2 described in eqn 3 of page1 of reference 1. Approximated as mentioned on
% page 3 of the reference. 
%
% Code by
% Aditya Kotikalpudi
% University of Minnesota
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I2 = zeros(size(u1));
I2_0 = zeros(size(u1));
I2_neg = zeros(size(u1));
u_temp1 = zeros(size(u1));
u1_temp2 = zeros(size(u1));
k1_temp1 = zeros(size(u1));
k1_temp2 = zeros(size(u1));

%% calculate I2
ind1 = find(u1>=0);
u_temp1(ind1) = u1(ind1);   % select elements in u1 > 0
k1_temp1(ind1) = k1(ind1);
I2_temp1 = getI2pos(u_temp1,k1_temp1);
I2(ind1) = I2_temp1(ind1);
j = sqrt(-1);
% Calulate integral I2(ref 1, page 1, eq 3) for u1<0
% Method taken from ref 3, page 90, eq 275
ind2 = find(u1<0);
u1_temp2(ind2) = u1(ind2);
k1_temp2(ind2) = k1(ind2);

I2_0temp = getI2pos(zeros(size(u1)),k1_temp2);
I2_0(ind2) = I2_0temp(ind2);
I2_negtemp = getI2pos(-u1_temp2,k1_temp2);
I2_neg(ind2) = I2_negtemp(ind2);
I2(ind2) = (2*real(I2_0(ind2))) - real(I2_neg(ind2)) + (j*imag(I2_neg(ind2))); 
return
     
         
          