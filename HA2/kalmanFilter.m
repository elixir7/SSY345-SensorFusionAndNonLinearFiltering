function [X, P, V] = kalmanFilter(Y, x_0, P_0, A, Q, H, R)
%KALMANFILTER Filters measurements sequence Y using a Kalman filter. 
%
%Input:
%   Y           [m x N] Measurement sequence
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   A           [n x n] State transition matrix
%   Q           [n x n] Process noise covariance
%   H           [m x n] Measurement model matrix
%   R           [m x m] Measurement noise covariance
%
%Output:
%   x           [n x N] Estimated state vector sequence
%   P           [n x n x N] Filter error convariance
%   V           [m x N] inovation

% Parameters
N = size(Y,2);

n = length(x_0);
m = size(Y,1);

% Data allocation
X = zeros(n,N);
P = zeros(n,n,N);
V = zeros(1,N);

% 1. Predict the next state
% 2. Update the prediction with measurement
% 3. Save the update for next iteration

for k = 1:N
    [x_0, P_0] = linearPrediction(x_0, P_0, A, Q);
    [x_0, P_0, V_0] = linearUpdate(x_0, P_0, Y(:,k), H, R);
    V(:, k) = V_0;
    X(:,k) = x_0;
    P(:,:,k) = P_0;
end
end

