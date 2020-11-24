function [A0,A1,Ap,Alag,Blag,Clag] = rational_approx(AIC,b_poles,omega)
% Least-squares approximation of the of the matrix qhg (frequency response)
% in order to reduce the number of aerodynamic states. The rational
% transfer function is in the form:
%       q(s) = a_0 + a_1*s + a_2*s^2 + sum(a_(i+2)*s/(s+b_i))
%
% Input arguments
%       qhg: Frequency response matrix of generalized aerodynamic influence
%       coefficients
%       b_poles: Vector of aerodynamic lags (poles of rational transfer
%       function).      b_poles = [b_1,b_2,...,b_i,...,b_n] 
%       omega: Frequency grid.
%                       omega = [w_1,...,w_m]
% Output arguments
%       A_0, A_1, A_2: Coefficient matrices for the rational transfer
%       function approximation.
%       A_p: 3D-Array containing the coefficient matrices associated with
%       the first order terms.
%       [A_lag, B_lag, C_lag] : State-space representation of the first
%       order terms.omega
% 
% Ref: K. L. Roger, �Airplane math modeling methods for active control
%      design,� Advisory Group for Aerospace Research and Development
%      CP-228, pp. 4�1, 1977.
%ome bac

%% Set-up
NumOmega = length(omega);       % Number of frequencies
NumPoles = length(b_poles);     % Number of lag poles

%% Define least square problem
% If there is one omega = 0, then we want to make it a fixed point and fit
% rest of the data with least square. Otherwise all the frequencies are
% used to setup the least square problem
idx_omege0 = find(omega==0);  % Find index of omega==0

% Error if more than one omega==0 found
if length(idx_omege0)>1
    error('Expecting at most one omega values of zero')
end

% Define temporary varables
denm = ((ones(NumOmega,1)*b_poles).^2+(omega'*ones(1,NumPoles)).^2);
numr = (omega'.^2)*ones(1,NumPoles);
numi = omega'*ones(1,NumPoles).*(ones(NumOmega,1)*b_poles);

% Define A matrix for least square problem x = A\b
if isempty(idx_omege0)      % If all omegas are nonzero
    A = [ones(NumOmega,1) zeros(NumOmega,1) numr./denm;
         zeros(NumOmega,1) omega' numi./denm];
else    % If one omega is zero
    a0 = AIC(:,:,idx_omege0);
    if ~isreal(a0)
        error('Expecting real values of AIC for omega=0') 
    end
    A = [zeros(NumOmega,1) numr./denm;
           omega' numi./denm];
end


% Obtain b and solve X = A\b
X = zeros([size(AIC,1),size(AIC,2),NumPoles+2]);

for irow = 1:size(AIC,1)
    for icol = 1:size(AIC,2)
        
        % Define a diagonal scaling factor to improve conditioning
        diag_scale =  diag(([1./abs(squeeze(AIC(irow,icol,:)));...
            1./abs(squeeze(AIC(irow,icol,:)))]));
        diag_scale(isinf(diag_scale)) = 1;
        A_scale = diag_scale*A;
        
        
        if isempty(idx_omege0)
            b_scale = diag_scale*[real(squeeze(AIC(irow,icol,:)));...
                imag(squeeze(AIC(irow,icol,:)))];
            X(irow,icol,:) = A_scale\b_scale;
        else
            b_scale = diag_scale*[real(squeeze(AIC(irow,icol,:))) ... 
                - a0(irow,icol,:);   ...
                imag(squeeze(AIC(irow,icol,:)))];
            X(irow,icol,:) = [a0(irow,icol,:); A_scale\b_scale];
        end
    end
end

A0 = X(:,:,1);
A1 = X(:,:,2);
Ap = X(:,:,3:end);

% %% Convert Lag transsfer functions to state space
% % Note that instead of the entire term A_p * (s)/(s+b), only A_p*1/(s+b) is
% % converted to state space. The simulink model takes care of the factor of
% % 's' in the numerator.
Alag = [];
Blag = [];
Clag = [];

for ii = 1:NumPoles
    Alag = blkdiag(Alag,-b_poles(ii)*eye(size(AIC,2)));
    Blag = [Blag; eye(size(AIC,2))];
    Clag = [Clag Ap(:,:,ii)];
end

end