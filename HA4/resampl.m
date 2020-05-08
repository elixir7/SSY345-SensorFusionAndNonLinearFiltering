function [Xr, Wr, j] = resampl(X, W)
%RESAMPLE Resample particles and output new particles and weights.
% resampled particles. 
%
%   if old particle vector is x, new particles x_new is computed as x(:,j)
%
% Input:
%   X   [n x N] Particles, each column is a particle.
%   W   [1 x N] Weights, corresponding to the samples
%
% Output:
%   Xr  [n x N] Resampled particles, each corresponding to some particle 
%               from old weights.
%   Wr  [1 x N] New weights for the resampled particles.
%   j   [1 x N] vector of indices refering to vector of old particles

    N = length(W);
    
    %1. Sample 
    Rs = rand(1, N);
    is = discretize(Rs, cumsum([0 W]));
    Xr = X(:, is);
    
    %2. Calculate weights
    Wr = 1/N * ones(1, N);
    
    %3. Index 
   % j = 2:(N+1);
    %j(end) = 1;
   % j = 1:N;
    j = is;
end