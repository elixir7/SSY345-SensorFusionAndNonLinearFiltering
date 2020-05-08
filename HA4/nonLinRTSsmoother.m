function [xs, Ps, xf, Pf, xp, Pp] = ...
    nonLinRTSsmoother(Y, x_0, P_0, f, T, Q, S, h, R, sigmaPoints, type)
% NONLINRTSSMOOTHER Filters measurement sequence Y using a 
% non-linear Kalman filter. 
%
%Input:
%   Y           [m x N] Measurement sequence for times 1,...,N
%   x_0         [n x 1] Prior mean for time 0
%   P_0         [n x n] Prior covariance
%   f                   Motion model function handle
%   T                   Sampling time
%   Q           [n x n] Process noise covariance
%   S           [n x N] Sensor position vector sequence
%   h                   Measurement model function handle
%   R           [n x n] Measurement noise covariance
%   sigmaPoints Handle to function that generates sigma points.
%   type        String that specifies type of non-linear filter/smoother
%
%Output:
%   xf          [n x N]     Filtered estimates for times 1,...,N
%   Pf          [n x n x N] Filter error convariance
%   xp          [n x N]     Predicted estimates for times 1,...,N
%   Pp          [n x n x N] Filter error convariance
%   xs          [n x N]     Smoothed estimates for times 1,...,N
%   Ps          [n x n x N] Smoothing error convariance

    N = size(Y,2);
    
    % Forward Kalman filter
    for k = 1:N
        % Prediction
        [x_0, P_0] = nonLinKFprediction(x_0, P_0, f, T, Q, sigmaPoints, type);
        xp(:,k) = x_0;
        Pp(:,:,k) = P_0;
        
        % Update
        [x_0, P_0] = nonLinKFupdate(x_0, P_0, Y(:,k), S(:,k), h, R, sigmaPoints, type);
        xf(:,k) = x_0;
        Pf(:,:,k) = P_0;
    end
    
    % RTS backward smoother
    xs(:,N) = xf(:,N); 
    Ps(:,:,N) = Pf(:,:,N);
    for k = (N-1):-1:1
        [xs(:,k) Ps(:,:,k)] = nonLinRTSSupdate(xs(:,k+1), Ps(:,:,k+1), xf(:,k), Pf(:,:,k), xp(:,k+1), Pp(:,:,k+1), f, T, sigmaPoints, type);    
    end

end