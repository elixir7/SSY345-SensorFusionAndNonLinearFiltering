function X = genNonLinearStateSequence(x_0, P_0, f, Q, N)
%GENLINEARSTATESEQUENCE generates an N+1-long sequence of states using a 
%    Gaussian prior and a linear Gaussian process model
%
%Input:
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   f           Motion model function handle
%               [fx,Fx]=f(x) 
%               Takes as input x (state), 
%               Returns fx and Fx, motion model and Jacobian evaluated at x
%               All other model parameters, such as sample time T,
%               must be included in the function
%   Q           [n x n] Process noise covariance
%   N           [1 x 1] Number of states to generate
%
%Output:
%   X           [n x N+1] State vector sequence
%

%Generate initial state from prior
x0 = mvnrnd(x_0, P_0, 1)';

% Generate state sequence by letting the initial state propagate
X = zeros(length(x_0), N+1);  % State sequence pre allocation
X(:,1) = x0;                  % Add initial state

for k = 2:N+1
    X(:,k) = mvnrnd( f(X(:,k-1)), Q, 1);
end


end