function Y = genNonLinearMeasurementSequence(X, h, R)
%GENNONLINEARMEASUREMENTSEQUENCE generates ovservations of the states 
% sequence X using a non-linear measurement model.
%
%Input:
%   X           [n x N+1] State vector sequence
%   h           Measurement model function handle
%               [hx,Hx]=h(x) 
%               Takes as input x (state) 
%               Returns hx and Hx, measurement model and Jacobian evaluated at x
%   R           [m x m] Measurement noise covariance
%
%Output:
%   Y           [m x N] Measurement sequence
%

% Pre allocate for measurement
Y = zeros(size(R,1), size(X,2)-1);

% Create measurement data from all states except the first
for k = 1:size(Y,2)
   Y(:,k) = mvnrnd(h(X(:,k+1)), R, 1); 
end

end