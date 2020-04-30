function [x, P] = nonLinKFupdate(x, P, y, h, R, type)
%NONLINKFUPDATE calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   y           [m x 1] measurement vector
%   h           Measurement model function handle
%               [hx,Hx]=h(x) 
%               Takes as input x (state), 
%               Returns hx and Hx, measurement model and Jacobian evaluated at x
%               Function must include all model parameters for the particular model, 
%               such as sensor position for some models.
%   R           [m x m] Measurement noise covariance
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] updated state mean
%   P           [n x n] updated state covariance
%
m = length(y);
n = length(x);
switch type
    case 'EKF'
       [hx Hx] = h(x);
       S = Hx*P*Hx' + R;
       K = P*Hx'*inv(S);
       x = x + K*(y - hx);
       P = P - K*S*K';
        
    case 'UKF'
       [SP, W] = sigmaPoints(x,P,type);
       for i = 1:length(SP)
           hx(:,i) = h(SP(:,i));
       end
       
       % Pre allocation
       y_pred = 0;
       P_xy = 0;
       S = R;
       for i = 1:(2*n+1)
           y_pred = y_pred + W(i).*hx(:,i);
       end
       for i = 1:(2*n+1)
           S = S +  (hx(:,i)-y_pred)*(hx(:,i)-y_pred)' .* W(i);
           P_xy = P_xy + ( (SP(:,i)-x) * (hx(:,i)-y_pred)' ) .*W(i);
       end
       
       %Measurement update
       x = x + P_xy*inv(S)*(y - y_pred)
       P = P - P_xy*inv(S)*P_xy';

        % Make sure the covariance matrix is semi-definite
        if min(eig(P))<=0
            [v,e] = eig(P, 'vector');
            e(e<0) = 1e-4;
            P = v*diag(e)/v;
        end
        
    case 'CKF'
       [SP, W] = sigmaPoints(x,P,type);
       for i = 1:length(SP)
           hx(:,i) = h(SP(:,i));
       end
       
       % Pre allocation
       y_pred = 0;
       P_xy = 0;
       S = R;
       for i = 1:(2*n)
           y_pred = y_pred + W(i).*hx(:,i);
       end
       for i = 1:(2*n)
           S = S +  (hx(:,i)-y_pred)*(hx(:,i)-y_pred)' .* W(i);
           P_xy = P_xy + ( (SP(:,i)-x) * (hx(:,i)-y_pred)' ) .*W(i);
       end
       
       %Measurement update
       x = x + P_xy*inv(S)*(y - y_pred)
       P = P - P_xy*inv(S)*P_xy';
        
    otherwise
        error('Incorrect type of non-linear Kalman filter')
end

end