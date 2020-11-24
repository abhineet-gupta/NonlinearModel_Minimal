function I1 = getI1(u1,k1)
%% Function to get I1
% ref 1: Albano and Rodden - A Doublet-Lattic Method for Calculating 
%        Lift Distributions on Oscillating Surfaces in Subsonic Flows
%
%ref 3: Blair, Max. A Compilation of the mathematics leading to the doublet 
%       lattice method. No. WL-TR-92-3028. WRIGHT LAB WRIGHT-PATTERSON AFB OH, 1992.
%
% I1 described in eqn 3 of page1 of reference 1. Approximated as shown in
% page 3 of the reference. 
%
% Code by
% Aditya Kotikalpudi
% University of Minnesota
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I1 = zeros(size(u1));
I1_0 = zeros(size(u1));
I1_neg = zeros(size(u1));
u_temp1 = zeros(size(u1));
u1_temp2 = zeros(size(u1));
k1_temp1 = zeros(size(u1));
k1_temp2 = zeros(size(u1));
%% evaluate I1 for u1>0
ind1 = find(u1>=0);
u_temp1(ind1) = u1(ind1);   % select elements in u1 > 0
k1_temp1(ind1) = k1(ind1);
I1_temp1 = getI1pos(u_temp1,k1_temp1);
I1(ind1) = I1_temp1(ind1);
j = sqrt(-1);

%% evaluate I1 for u1<0
% Method taken from ref 3, page 90, eq 275
ind2 = find(u1<0);
u1_temp2(ind2) = u1(ind2);
k1_temp2(ind2) = k1(ind2);

I1_0temp = getI1pos(zeros(size(u1)),k1_temp2);
I1_0(ind2) = I1_0temp(ind2);
I1_negtemp = getI1pos(-u1_temp2,k1_temp2);
I1_neg(ind2) = I1_negtemp(ind2);
I1(ind2) = (2*real(I1_0(ind2))) - real(I1_neg(ind2)) + (j*imag(I1_neg(ind2)));

return
