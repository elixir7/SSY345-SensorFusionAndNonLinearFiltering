function [xs, Ps] = nonLinRTSSupdate(xs_kplus1, ...
                                     Ps_kplus1, ...
                                     xf_k, ... 
                                     Pf_k, ...
                                     xp_kplus1, ...
                                     Pp_kplus1, ...
                                     f, ...
                                     T, ...
                                     sigmaPoints, ...
                                     type)
%NONLINRTSSUPDATE Calculates mean and covariance of smoothed state
% density, using a non-linear Gaussian model.
%
%Input:
%   xs_kplus1   Smooting estimate for state at time k+1
%   Ps_kplus1   Smoothing error covariance for state at time k+1
%   xf_k        Filter estimate for state at time k
%   Pf_k        Filter error covariance for state at time k
%   xp_kplus1   Prediction estimate for state at time k+1
%   Pp_kplus1   Prediction error covariance for state at time k+1
%   f           Motion model function handle
%   T           Sampling time
%   sigmaPoints Handle to function that generates sigma points.
%   type        String that specifies type of non-linear filter/smoother
%
%Output:
%   xs          Smoothed estimate of state at time k
%   Ps          Smoothed error convariance for state at time k

    P_k_kplus1 = 0;
    n = length(xf_k);

    if(type == "EKF")
        [fx, Fx] = f(xf_k, T);
        P_k_kplus1 = Pf_k * Fx';
    elseif(type == "UKF")
        [SP,W] = sigmaPoints(xf_k, Pf_k, 'UKF');
        %Transform sigma points
        for i = 1:length(SP)
            fx(:,i) = f(SP(:,i), T);
        end

        P = zeros(n,n);
        for i = 1:(2*n+1)
            P = P + ( (SP(:,i)-xf_k) * (fx(:,i)-xp_kplus1)' ) .*W(i);
        end
        P_k_kplus1 = P;
    elseif(type == "CKF")
        [SP,W] = sigmaPoints(xf_k, Pf_k, 'CKF');
        %Transform sigma points
        for i = 1:length(SP)
            fx(:,i) = f(SP(:,i), T);
        end

        P = zeros(n,n);
        for i = 1:(2*n)
            P = P + ( (SP(:,i)-xf_k) * (fx(:,i)-xp_kplus1)' ) .*W(i);
        end
        P_k_kplus1 = P;
    else
        error("type is not of EKF,UKF or CKF")
    end

    G_k = P_k_kplus1*inv(Pp_kplus1);
    xs = xf_k + G_k*(xs_kplus1 - xp_kplus1);
    Ps = Pf_k - G_k*(Pp_kplus1 - Ps_kplus1)*G_k';

end