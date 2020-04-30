function [fx, Fx] = coordinatedTurnMotion(x, T)
%COORDINATEDTURNMOTION calculates the predicted state using a coordinated
%turn motion model, and also calculated the motion model Jacobian
%
%Input:
%   x           [5 x 1] state vector
%   T           [1 x 1] Sampling time
%
%Output:
%   fx          [5 x 1] motion model evaluated at state x
%   Fx          [5 x 5] motion model Jacobian evaluated at state x
%
% NOTE: the motion model assumes that the state vector x consist of the
% following states:
%   px          X-position
%   py          Y-position
%   v           velocity
%   theta       heading
%   omega       turn-rate

% For easy visualisation, pull out the members and assign appropriate names
px = x(1);
py = x(2);
v = x(3);
theta = x(4);
omega = x(5);


% Next time step x(k) = f(x(k-1)) + q(k-1)
fx = [
  px + T*v*cos(theta);
  py + T*v*sin(theta);
  v;
  theta + T*omega;
  omega;
];

%Check if the Jacobian is requested by the calling function
if nargout > 1
    Fx = [
        1   0   T*cos(theta)    -T*v*sin(theta) 0;
        0   1   T*sin(theta)    T*v*cos(theta)  0;
        0   0   1               0               0;
        0   0   0               1               T;
        0   0   0               0               1
    ];
end

end