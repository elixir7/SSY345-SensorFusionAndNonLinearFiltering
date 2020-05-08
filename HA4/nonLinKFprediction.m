function [x, P] = nonLinKFprediction(x, P, f, T, Q, sigmaPoints, type)
%NONLINKFPREDICTION calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   f           Motion model function handle
%               [fx,Fx]=f(x,T) 
%               Takes as input x (state) and T sample time 
%               Returns fx and Fx, motion model and Jacobian evaluated at x
%   Q           [n x n] Process noise covariance
%   sigmaPoints Handle to function that generates sigma points.
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] predicted state mean
%   P           [n x n] predicted state covariance
%

% Prediction
%[fx, Fx] = f(x);
n = size(x,1);
switch type
    case 'EKF'
        
        % 1. Predict state using non linear function
        % 2. Estimate the gaussian covariance using a linearization
        % Note: The acctual distrubution is not gaussian but we approximate 
        %       to a gaussian to simplify.
        [fx, Fx] = f(x,T);
        x = fx;
        P = Fx*P*Fx'+ Q;
        
    case 'UKF'
        % 1. Get sigma points
        % 2. Transform SP using non linear transform
        % 3. Calculate estimated mean using transformed SP
        % 4. Calculate estimated covar using transformed SP and mean
        
        [SP,W] = sigmaPoints(x, P, 'UKF');
        %Transform sigma points
        for i = 1:length(SP)
            fx(:,i) = f(SP(:,i), T);
        end
        
        x = 0;
        for i = 1:(2*n+1)
            x = x + W(i).*fx(:,i);
        end
        
        P = zeros(n,n);
        for i = 1:(2*n+1)
            P = P + ( (fx(:,i)-x) * (fx(:,i)-x)' ) .*W(i);
        end
        P = P + Q;
        
        % Make sure the covariance matrix is semi-definite
        if min(eig(P))<=0
            [v,e] = eig(P, 'vector');
            e(e<0) = 1e-4;
            P = v*diag(e)/v;
        end
            
    case 'CKF'
        % 1. Get sigma points
        % 2. Transform SP using non linear transform
        % 3. Calculate estimated mean using transformed SP
        % 4. Calculate estimated covar using transformed SP and mean
        
        [SP,W] = sigmaPoints(x, P, 'CKF');
        %Transform sigma points
        for i = 1:length(SP)
            fx(:,i) = f(SP(:,i), T);
        end
        
        x = 0;
        for i = 1:(2*n)
            x = x + W(i).*fx(:,i);
        end
        
        P = zeros(n,n);
        for i = 1:(2*n)
            P = P + W(i).*((fx(:,i)-x) * (fx(:,i)-x)');
        end
        P = P + Q;
        
    otherwise
        error('Incorrect type of non-linear Kalman filter')
end

end