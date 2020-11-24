function I2pos = getI2pos(u1,k1)
% this function gets I2 integral for non-planar body solutions, for u1>0
% I2 = I2_1 + I2_2
% Expressions for I2_1 & I2_2 have been derived using the same
% approximations as those for I1 
a1 = 0.101;
a2 = 0.899;
a3 = 0.09480933;
b1 = 0.329;
b2 = 1.4067;
b3 = 2.90;
i = sqrt(-1);
eiku = exp(-i.*k1.*u1);


I2_1 = getI1pos(u1,k1);
%I2_2_1 = 0;
I2_2_1 = a1.*exp(-(b1+i*k1).*u1)./((b1+i*k1).^2) + a2.*exp(-(b2+i*k1).*u1)./((b2+i*k1).^2) ...
         + ((a3.*exp(-(b3+i*k1).*u1)./(((b3+i*k1).^2 + pi^2).^2)).*(pi*((pi*sin(pi*u1)) - ((b3+i*k1).*cos(pi*u1))) ...
         - ((b3+i*k1).*(pi*cos(pi*u1) + ((b3+i*k1).*sin(pi*u1))))));  
I2_2 = (eiku.*(u1.^3)./((1+u1.^2).^(3/2)) - getI1pos(u1,k1) - ...
       eiku.*u1./sqrt((1+u1.^2)))/3 - (k1.*k1.*I2_2_1/3);
  
I2pos = I2_1 + I2_2;
return