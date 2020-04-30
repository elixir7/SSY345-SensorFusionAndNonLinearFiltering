function [xf, Pf, xp, Pp] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q, h, R, type)
    % NONLINEARKALMANFILTER Filters measurement sequence Y using a 
    % non-linear Kalman filter. 
    %
    % Input:
    %   Y           [m x N] Measurement sequence for times 1,...,N
    %   x_0         [n x 1] Prior mean for time 0
    %   P_0         [n x n] Prior covariance
    %   f                   Motion model function handle
    %                       [fx,Fx]=f(x) 
    %                       Takes as input x (state) 
    %                       Returns fx and Fx, motion model and Jacobian evaluated at x
    %   Q           [n x n] Process noise covariance
    %   h                   Measurement model function handle
    %                       [hx,Hx]=h(x,T) 
    %                       Takes as input x (state), 
    %                       Returns hx and Hx, measurement model and Jacobian evaluated at x
    %   R           [m x m] Measurement noise covariance
    %
    % Output:
    %   xf          [n x N]     Filtered estimates for times 1,...,N
    %   Pf          [n x n x N] Filter error convariance
    %   xp          [n x N]     Predicted estimates for times 1,...,N
    %   Pp          [n x n x N] Filter error convariance

    % Parameters
    N = size(Y, 2);

    n = length(x_0);
    m = size(Y, 1);

    % Data allocation
    xp = zeros(n, N);
    Pp = zeros(n, n, N);
    xf = zeros(n, N);
    Pf = zeros(n, n, N);

    for i = 1:N
        % x_0 and P_0 holds last states values, i.e.
        % x(k-1) and P(k-1). With this we make a
        % prediction for x(k | k-1).
        [x_0, P_0] = nonLinKFprediction(x_0, P_0, f, Q, type);
        
        % Store x(k | k-1) and P(k | k-1) into array.
        xp(:, i) = x_0;
        Pp(:, :, i) = P_0;
        
        % Now, x_0 and P_0 holds x(k) and P(k) (the prediciton).
        % We will now preform an update, i.e. x(k | k).
        [x_0, P_0] = nonLinKFupdate(x_0, P_0, Y(:, i), h, R, type);
        
        % Store x(k | k) and P(k | k) into array.
        xf(:, i) = x_0;
        Pf(:, :, i) = P_0;
    end

end
