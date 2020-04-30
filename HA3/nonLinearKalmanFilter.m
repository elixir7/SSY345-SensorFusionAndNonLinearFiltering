function [xf, Pf, xp, Pp] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q, h, R, type)
%NONLINEARKALMANFILTER Filters measurement sequence Y using a 
% non-linear Kalman filter. 
%
%Input:
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
%Output:
%   xf          [n x N]     Filtered estimates for times 1,...,N
%   Pf          [n x n x N] Filter error convariance
%   xp          [n x N]     Predicted estimates for times 1,...,N
%   Pp          [n x n x N] Filter error convariance
%

% Parameters
N = size(Y,2);

%n = length(x_0);
%m = size(Y,1);

% Data allocation
%X = zeros(n,N);
%P = zeros(n,n,N);
%V = zeros(1,N);

% 1. Predict the next state
% 2. Update the prediction with measurement
for k = 1:N
    % Prediction
    [x_0, P_0] = nonLinKFprediction(x_0, P_0, f, Q, type);
    xp(:,k) = x_0;
    Pp(:,:,k) = P_0;
    
    % Update
    [x_0, P_0] = nonLinKFupdate(x_0, P_0, Y(:,k), h, R, type);
    xf(:,k) = x_0;
    Pf(:,:,k) = P_0;
    
    
    
end


end