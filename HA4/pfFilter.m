function [xfp, Pfp, Xp, Wp] = pfFilter(x_0, P_0, Y, proc_f, proc_Q, meas_h, meas_R, ...
                             N, bResample, plotFunc)
%PFFILTER Filters measurements Y using the SIS or SIR algorithms and a
% state-space model.
%
% Input:
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   Y           [m x K] Measurement sequence to be filtered
%   proc_f      Handle for process function f(x_k-1)
%   proc_Q      [n x n] process noise covariance
%   meas_h      Handle for measurement model function h(x_k)
%   meas_R      [m x m] measurement noise covariance
%   N           Number of particles
%   bResample   boolean false - no resampling, true - resampling
%   plotFunc    Handle for plot function that is called when a filter
%               recursion has finished.
% Output:
%   xfp         [n x K] Posterior means of particle filter
%   Pfp         [n x n x K] Posterior error covariances of particle filter
%   Xp          [n x N x K] Particles for posterior state distribution in times 1:K
%   Wp          [N x K] Non-resampled weights for posterior state x in times 1:K


    % If you want to be a bit fancy, then only store and output the particles if the function
    % is called with more than 2 output arguments.
    K = size(Y,2);
    n = length(x_0);

    % Allocate memory
    xfp = zeros(n, K);
    Pfp = zeros(n,n,K);
    Xp = zeros(n,N,K);
    Wp = zeros(N, K);


    %1. Create particles
    %2. For each time instant
        %2.1 Filter each particle
        %2.2 Resample if asked
        %2.3 Compute state density with montecarlo (with weights)
        %2.4 Approximate mean and covar.


    X0 = mvnrnd(x_0, P_0, N)';
    W0 = 1/N * ones(1,N);

    j = 1:N;
    
    for k = 1:K
        Xmin1 = X0;
        
        [X0, W0] = pfFilterStep(X0, W0, Y(:,k), proc_f, proc_Q, meas_h, meas_R);
        
        % Plot
        if nargin >= 10 
            plotFunc(k, X0, Xmin1, W0, j)
        end
        
        % Resample
        if bResample
            [X0, W0, j] = resampl(X0, W0);
        end
        Xp(:,:,k) = X0;
        Wp(:, k) = W0;

        % Mean
        xfp(:,k) = sum(W0.*X0, 2);

        % Covar
        Pfp(:,:,k) = (W0.*(X0 - xfp(:,k)))*(X0 - xfp(:,k))';
        
        
    end

end