function [estimatedPos, tdoa] = solveTDOA(toa, gNBPos, x0)
    c = physconst('LightSpeed');
    % TDOA computation
    tdoa = toa(2:end)-toa(1);

    % ANALYTICAL METHOD
    % Plot gNB and UE positions
    % Selects reference base station (first one)
    S1        = gNBPos(1,:);           % 1 x 3 reference base station position
    Sk        = gNBPos(2:end,:);       % (N-1) x 3 remaining base station positions
    deltaPos  = Sk - S1;               % (N-1) x 3 relative position vectors
    dkt       = c * tdoa(:);           % (N-1) x 1 range differences
    
    % Build matrix A and vector b
    A = [ deltaPos, -dkt ];  
    b = ( sum(Sk.^2, 2) ...            % Row-wise squared norm of Sk
        - sum(S1.^2) ...               % Squared norm of S1
        - dkt.^2 ) / 2;                % Subtract dk squared and divide by 2
        
    % Minimum-norm least squares solver
    % Solved using lsqlin with constraints
    f  = [];                           % No additional quadratic cost function
    lb = [-Inf, -Inf, -5, -Inf];       % Lower bounds
    ub = [ Inf,  Inf,  5,  Inf];       % Upper bounds
    
    % Solve the constrained system
    % x0 = [mean(UEPos), 0];           % Initial guess: [x y z r1]
                                       % Mean BS position and initial r1 = 0
    % x0 = UEPos;                      % Alternative initial guess
    options = optimoptions('lsqlin', 'Algorithm', 'interior-point', 'Display', 'off');
    x_est = lsqlin(A, b, [], [], [], [], lb, ub, x0, options);

    % Extract only the estimated position [x, y, z]
    estimatedPos = x_est(1:3);
end
