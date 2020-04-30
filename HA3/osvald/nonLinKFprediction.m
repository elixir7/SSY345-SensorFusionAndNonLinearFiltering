function [x, P] = nonLinKFprediction(x_0, P_0, f, Q, type)
    % NONLINKFPREDICTION calculates mean and covariance of predicted state
    %   density using a non-linear Gaussian model.
    %
    % Input:
    %   x_0         [n x 1] Prior mean
    %   P_0         [n x n] Prior covariance
    %   f           Motion model function handle
    %               [fx,Fx]=f(x) 
    %               Takes as input x (state), 
    %               Returns fx and Fx, motion model and Jacobian evaluated at x
    %               All other model parameters, such as sample time T,
    %               must be included in the function
    %   Q           [n x n] Process noise covariance
    %   type        String that specifies the type of non-linear filter
    %
    % Output:
    %   x           [n x 1] predicted state mean
    %   P           [n x n] predicted state covariance
    %

    if (type == 'EKF')

        [fx, Fx] = f(x_0);
        x = fx;
        P = Fx * P_0 * Fx' + Q;
            
    elseif (type == 'UKF' | type == 'CKF')
        % UKF and CKF is the same method.
        % Only the sigma points differ.
        
        n = length(x_0);
        
        [SP, W] = sigmaPoints(x_0, P_0, type);
        N = length(SP);
        
        % Allocate memory and calculate
        % the non-linear transform.
        fSP = zeros(n, N);
        for i = 1:N
            fSP(:, i) = f(SP(:, i));
        end
        
        % Make the prediction
        x = sum(W' .* fSP')';
        e = fSP - x;
        P = (W .* e) * e' + Q;

        % Make sure the covariance matrix is semi-definite
        % (only a concern for UKF)
        if min(eig(P))<=0
            [v,e] = eig(P, 'vector');
            e(e<0) = 1e-4;
            P = v*diag(e)/v;
        end
            
    else
        error('Incorrect type of non-linear Kalman filter');
    end

end
