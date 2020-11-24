function StateVals = SetInitial(Vinf,AircraftType)
%% SetTrimTarget     Sets trim targets to be used as trim initial condition

% switch FlightCond
% case 1
    % Set initial control surface positions
%     StateVals.CSPos = pi/180*[-0.6, -0.6, -0.6, -0.6, ...
%                             -0.6, -0.6, -0.6, -0.6];
    StateVals.CSPos = [0, 0, 0.1, 0, ...
                            0, 0, 0.1, 0];
    
    % Set initial control surface velocities
    StateVals.CSVel = [0, 0, 0, 0, 0, 0, 0, 0];
    
    % Set initial Eular angles
    StateVals.Eular = pi/180*[0, 2.5, 0];
    
    % Set initial angular rates
    StateVals.AngRates = pi/180*[0, 0, 0];
    
    % Set initial body frame velocities
    StateVals.Vb = [Vinf 0 1];
    
    % Set initial position in earth frame
    % (x>Towards nose, y>Towards starboard, z> Towards earth) 
    StateVals.Xe = [0 0 -100];
    
    % Set initial structural deflections
    StateVals.Eta = [ 0 0 0 0 0];
    
    % Set initial structural deflections
    StateVals.EtaDot = [ 0 0 0 0 0];
    
    % Set initail lage states
    if strcmp(AircraftType,'DLM')
        StateVals.Lag = zeros(1,888*2);
    end
% end
