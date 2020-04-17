function Y = genLinearMeasurementSequence(X, H, R)
%GENLINEARMEASUREMENTSEQUENCE generates a sequence of observations of the state 
% sequence X using a linear measurement model. Measurement noise is assumed to be 
% zero mean and Gaussian.
%
%Input:
%   X           [n x N+1] State vector sequence. The k:th state vector is X(:,k+1)
%   H           [m x n] Measurement matrix
%   R           [m x m] Measurement noise covariance
%
%Output:
%   Y           [m x N] Measurement sequence
%

% your code here
% Create measurement vector
Y = zeros(size(H,1), size(X,2)-1);

% Create measurement data from all states except the first
for k = 1:size(Y,2)
   Y(:,k) = mvnrnd(H*X(:,k+1), R, 1); 
end


end