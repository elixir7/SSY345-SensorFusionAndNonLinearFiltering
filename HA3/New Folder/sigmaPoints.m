function [SP,W] = sigmaPoints(x, P, type)
% SIGMAPOINTS computes sigma points, either using unscented transform or
% using cubature.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%
%Output:
%   SP          [n x 2n+1] UKF, [n x 2n] CKF. Matrix with sigma points
%   W           [1 x 2n+1] UKF, [1 x 2n] UKF. Vector with sigma point weights 
%
n = size(x,1);
P_sqrt = sqrtm(P);

switch type        
    case 'UKF'
        num_points = 2*n + 1;
        
        SP = zeros(n, num_points);
        SP(:,1) = x;
        
        W = zeros(1, num_points);
        W0 = 1 - n/3;
        
        W(1) = W0;
        SP(:,1) = x;
        
        for i = 1:n
           W(i+1) = (1-W0) / (2*n); 
           W(i+1+n) = (1-W0) / (2*n); 
           
           SP(:,i+1)   = x + sqrt(n/(1-W0))*P_sqrt(:,i);
           SP(:,i+1+n) = x - sqrt(n/(1-W0))*P_sqrt(:,i);
        end
        
            
    case 'CKF'
        num_points = 2*n;
        
        SP = zeros(n, num_points);
        W = zeros(1, num_points);
        
        for i = 1:n
            W(i)      = 1/(2*n);
            W(i+n)    = W(i); 
            SP(:,i)   = x + sqrt(n)*P_sqrt(:,i);
            SP(:,i+n) = x - sqrt(n)*P_sqrt(:,i);
        end
        
    otherwise
        error('Incorrect type of sigma point')
end

end