function [x, P] = nonLinKFupdate(x_0, P_0, y_0, h, R, type)
    % NONLINKFUPDATE calculates mean and covariance of predicted state
    %   density using a non-linear Gaussian model.
    %
    % Input:
    %   x_0         [n x 1] Prior mean
    %   P_0         [n x n] Prior covariance
    %   y_0         [m x 1] measurement vector
    %   h           Measurement model function handle
    %               [hx,Hx]=h(x) 
    %               Takes as input x (state), 
    %               Returns hx and Hx, measurement model and Jacobian evaluated at x
    %               Function must include all model parameters for the particular model, 
    %               such as sensor position for some models.
    %   R           [m x m] Measurement noise covariance
    %   type        String that specifies the type of non-linear filter
    %
    % Output:
    %   x           [n x 1] updated state mean
    %   P           [n x n] updated state covariance

    if (type == 'EKF')
            
        % Deriviation of extended kalman filter
        % update operation comes from the L6.
        [hx, Hx] = h(x_0);

        S = Hx * P_0 * Hx' + R;
        K = P_0 * Hx' * inv(S);
        v = K * (y_0 - hx);

        x = x_0 + v;
        P = P_0 - K * S * K';

    elseif (type == 'UKF' | type == 'CKF')
        % UKF and CKF is the same method.
        % Only the sigma points differ.

        m = length(y_0);
        
        [SP, W] = sigmaPoints(x_0, P_0, type);
        N = length(SP);
        
        % Allocate memory and calculate
        % the non-linear transform.
        hSP = zeros(m, N);
        for i = 1:N
            hSP(:, i) = h(SP(:, i));
        end
        
        % Make the update
        % y is y_hat
        y = sum(W' .* hSP')';
        e = hSP - y;
        Pxy = (W .* (SP - x_0)) * e';
        S = (W .* e) * e' + R;

        x = x_0 + Pxy * inv(S) * (y_0 - y);
        P = P_0 - Pxy * inv(S) * Pxy';

        % Make sure the covariance matrix is semi-definite
        % (only a concern for UKF)
        if min(eig(P))<=0
            [v,e] = eig(P, 'vector');
            e(e<0) = 1e-4;
            P = v*diag(e)/v;
        end
        
    else
        error('Incorrect type of non-linear Kalman filter')
    end

end
