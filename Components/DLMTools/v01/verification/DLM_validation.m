%% DLM code verification
% ref 1: Rodden,W., Taylor, P.F., McIntosh,S.C., "Improvements to the Doublet Lattice Method
%        in MSC/NASTRAN"
%
% code to verify the DLM code for a basic benchmark case where a single
% panel is defined and the corresponding kernel function distribution is
% plotted. The test case is based on ref 1, page no. 7.
% resulting plots need to be compared against Fig. 4.a and 4.b of ref. 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

pathname = pwd;
addpath([pathname(1:end-13) '\DLM_code']);

%% create a single panel, chord = 0.2, span = 1
chord = 0.2;
span = 1;
e = 0.5*span;
N = 1;

x01 = chord*0.5;
x02 = chord*0.5;
x03 = chord*0.5;

y01 = span*0.5;
y02 = 0;
y03 = -span*0.5;

z01 = 0;
z02 = 0;
z03 = 0;

cosGamma = 1;
sinGamma = 0;

k = 2*2/7;  % chord is to be taken as 7, refer to ref. 1
Mach = 0.8;

%% get Ki,Km,K0
Ki_w = getKappa(x01,y01,z01,cosGamma,sinGamma,k,Mach);
Ki_0 = getKappa(x01,y01,z01,cosGamma,sinGamma,0,Mach);
Ki = Ki_w - Ki_0;       
                        
% numerator of the singular kernel for this doublet
Km_w = getKappa(x02,y02,z02,cosGamma,sinGamma,k,Mach);
Km_0 = getKappa(x02,y02,z02,cosGamma,sinGamma,0,Mach);
Km = Km_w - Km_0;

% numerator of the singular kernel for this doublet
K0_w = getKappa(x03,y03,z03,cosGamma,sinGamma,k,Mach);
K0_0 = getKappa(x03,y03,z03,cosGamma,sinGamma,0,Mach);
K0 = K0_w - K0_0;


%% get A,B,C

e1 = abs(repmat(e,N,1));
A = (Ki-2*Km+K0)./(2*e1.^2);
B = (K0-Ki)./(2*e1);
C = Km;

%save mydata

%% plot the parabola
eta = -0.5:0.01:0.5;
eta2 = eta.^2;

K_total = A*eta2 + B*eta + C;

figure(1)
plot(eta,real(K_total),'b-o')
title('real part of Kernel function')

figure(2)
plot(eta,imag(K_total),'b-o')
title('Imaginary part of Kernel function')