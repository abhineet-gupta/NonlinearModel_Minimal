function D = getDLM(P0,P1,P2,P3,e,cav,k,Mach)
%
% This function solves equation 6 using equation 7 from ref 1. Implements 
% the parabolic distribution assumption for the numerator of the
% incremental Kernel function. 
% 
%
% ref 1: Albano and Rodden - A Doublet-Lattic Method for Calculating 
%        Lift Distributions on Oscillating Surfaces in Subsonic Flows
%
% ref 2: Rodden,Taylor and McIntosh, 'Further refinement of Doublet lattice
%        method
% chord refers to chord of the aero panel
% span refers to the span of the aero panel
% 
% Input
%
% P0 = downwash recieving location x-y pair (1/2span, 3/4chord)
% P1 = root doublet location x-y pair (0 span, 1/4chord)
% P2 = semi-span doublet location x-y pair (1/2span, 1/4chord)
% P3 = tip doublet location x-y pair (1span, 1/4chord)
% e = half span length of the aero panel
% cav = centerline chord of the aero panel
% k = k1 term from ref 1.: omega/U
%   omega = frequency of oscillation
%   U = freestream velocity
% M = Mach number
%
% Output
%
% D = normalwash factor for the aero panel in question due to the doublets
% at the locations defined by P1, P2 and P3
%
% Based on
% Original code by Frank R. Chavez, Iowa State University (~2002)
% Modified version by:
% Brian P. Danowsky
% P. Chase Schulze
% (c) Systems Technology, Inc. 2014
%
% Latest modification by 
% Aditya Kotikalpudi
% University of Minnesota
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(P0,1);

%% Get relevant inputs for calculating Kernel function 
%(Page 1, reference 1,beginning with eq. 2)

% root doublet
% define the x,y,z distance from this doublet to the recieving point at
% 1/2span 3/4chord location
x01= bsxfun(@minus,P0(:,1),P1(:,1)');
y01= bsxfun(@minus,P0(:,2),P1(:,2)');
z01= bsxfun(@minus,P0(:,3),P1(:,3)');

% semi-span doublet
% define the x,y,z distance from this doublet to the recieving point at
% 1/2span 3/4chord location
x02= bsxfun(@minus,P0(:,1),P2(:,1)');
y02= bsxfun(@minus,P0(:,2),P2(:,2)');
z02= bsxfun(@minus,P0(:,3),P2(:,3)');

% tip doublet
% define the x,y,z distance from this doublet to the recieving point at
% 1/2span 3/4chord location
x03= bsxfun(@minus,P0(:,1),P3(:,1)');
y03= bsxfun(@minus,P0(:,2),P3(:,2)');
z03= bsxfun(@minus,P0(:,3),P3(:,3)');

%% cos and sin angles of all doublet line dihedrals

% CosGamma = y/sqrt(z^2 + y^2)
cosGamma = (P3(:,2) - P1(:,2))./ ...
            (sqrt((P3(:,3)-P1(:,3)).^2 + (P3(:,2)-P1(:,2)).^2));
            
% SinGamma = z/sqrt(z^2 + y^2)
sinGamma = (P3(:,3) - P1(:,3))./(sqrt((P3(:,3)-P1(:,3)).^2 + ...
            (P3(:,2)-P1(:,2)).^2));

%% Kernel function (K) calculation
% Kappa (defined on page 3, reference 1) is calculated. The steady part of
% Kappa (i.e. reduced frequency = 0) is subtracted out and later
% compensated for by adding downwash effects from a VLM code. This ensures 
% that the doublet lattice code converges to VLM results under steady
% conditions. (Ref 2, page 3, equation 9)

% T1 calculation for BFF
%T1 = [ones(Nw,Nw) zeros(Nw,Nwl); zeros(Nwl,Nw) ones(Nwl,Nwl)];
% numerator of the singular kernel for this doublet
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

%% Parabolic approximation of incremental Kernel function (ref 1, equation 7)
% define terms used in the parabolic approximation

e1 = abs(repmat(e',N,1));
A = (Ki-2*Km+K0)./(2*e1.^2);
B = (K0-Ki)./(2*e1);
C = Km;

% define r1,n0,zeta0

cosGamma = repmat(cosGamma',N,1);
sinGamma = repmat(sinGamma',N,1);

n0 = (y02.*cosGamma) + (z02.*sinGamma);
zeta0 = -(y02.*sinGamma) + (z02.*cosGamma);

r2 = sqrt((n0.^2) + (zeta0.^2));



% normalwash matrix factor I

I = (A.*(2*e1))+((0.5*B+n0.*A).*log((r2.^2 - 2*n0.*e1 + e1.^2)./(r2.^2 + 2*n0.*e1 + e1.^2)))...
    +(((n0.^2 - zeta0.^2).*A+n0.*B+C)./abs(zeta0).*atan(2*e1.*abs(zeta0)./(r2.^2 - e1.^2)));

% limit when zeta -> 0
ind = find(zeta0==0);
I2 = ((A.*(2*e1))+((0.5*B+n0.*A).*log(((n0-e1)./(n0+e1)).^2))...
    +((n0.^2).*A+n0.*B+C).*((2*e1)./(n0.^2 - e1.^2)));

I(ind) = I2(ind);

% normalwash matrix
D = repmat(cav',N,1).*I/(pi*8);


return
