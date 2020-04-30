%% HA3 Analysys
%% 1 Approximations of mean and covariance
%% a
clear; clc; close all;
format short g;

s1.x0 = [120; 120];
s1.P = diag([5^2, 10^2]);

s2.x0 = [120; -20];
s2.P = diag([5^2, 10^2]);

s1_pos = [0; 100];
s2_pos = [100; 0];

sigma_phi = (0.1*pi/180);
R = diag([sigma_phi^2, sigma_phi^2]);


% Sample from the distrubution
N = 10000;
s1.x = mvnrnd(s1.x0, s1.P, N)';
s2.x = mvnrnd(s2.x0, s2.P, N)';
for i = 1:N
    s1.y(:,i) = dualBearingMeasurement(s1.x(:,i), s1_pos, s2_pos);
    s1.y(:,i) = s1.y(:,i) + mvnrnd([0;0], R, 1)';  
    
    s2.y(:,i) = dualBearingMeasurement(s2.x(:,i), s1_pos, s2_pos);
    s2.y(:,i) = s2.y(:,i) + mvnrnd([0;0], R, 1)';  
end

% Monte carlo mean & covar approximation
s1.y_mc = sum(s1.y, 2) / N;
s2.y_mc = sum(s2.y, 2)/ N;

s1.y_covar = 0;
s2.y_covar = 0;
for i = 1:N
    s1.y_covar = s1.y_covar + (s1.y(:,i)-s1.y_mc)*(s1.y(:,i)-s1.y_mc)';
    s2.y_covar = s2.y_covar + (s2.y(:,i)-s2.y_mc)*(s2.y(:,i)-s2.y_mc)';
end
s1.y_covar = s1.y_covar/N;
s2.y_covar = s2.y_covar / N;

%% b
h = @(x) dualBearingMeasurement(x, s1_pos, s2_pos);
level = 3; 
npoints = 5000;

% State density 1
% Mean
[s1.x_ekf, s1.P_ekf] = nonLinKFprediction(s1.x0, s1.P, h, R, 'EKF');
[s1.x_ukf, s1.P_ukf] = nonLinKFprediction(s1.x0, s1.P, h, R, 'UKF');
[s1.x_ckf, s1.P_ckf] = nonLinKFprediction(s1.x0, s1.P, h, R, 'CKF');
% Covariance
[s1.xy_mc]  = sigmaEllipse2D( s1.y_mc, s1.y_covar, level, npoints );
[s1.xy_ekf] = sigmaEllipse2D( s1.x_ekf, s1.P_ekf, level, npoints );
[s1.xy_ukf] = sigmaEllipse2D( s1.x_ukf, s1.P_ekf, level, npoints );
[s1.xy_ckf] = sigmaEllipse2D( s1.x_ckf, s1.P_ekf, level, npoints );


% State density 2
% Mean
[s2.x_ekf, s2.P_ekf] = nonLinKFprediction(s2.x0, s2.P, h, R, 'EKF');
[s2.x_ukf, s2.P_ukf] = nonLinKFprediction(s2.x0, s2.P, h, R, 'UKF');
[s2.x_ckf, s2.P_ckf] = nonLinKFprediction(s2.x0, s2.P, h, R, 'CKF');
% Covariance
[s2.xy_mc]  = sigmaEllipse2D( s2.y_mc, s2.y_covar, level, npoints );
[s2.xy_ekf] = sigmaEllipse2D( s2.x_ekf, s2.P_ekf, level, npoints );
[s2.xy_ukf] = sigmaEllipse2D( s2.x_ukf, s2.P_ukf, level, npoints );
[s2.xy_ckf] = sigmaEllipse2D( s2.x_ckf, s2.P_ckf, level, npoints );


% Plot scenario 1
figure(1); clf; hold on; grid on;
%Samples
plot(s1.y(1,:), s1.y(2,:), '.')
options = {'MarkerSize', 10, 'LineWidth', 1};
% Mean
plot(s1.y_mc(1), s1.y_mc(2), 'yo', options{:});
plot(s1.x_ekf(1), s1.x_ekf(2), 'mo', options{:});
plot(s1.x_ukf(1), s1.x_ukf(2), 'go', options{:});
plot(s1.x_ckf(1), s1.x_ckf(2), 'ro', options{:});
% 3-Sigma Ellipse
plot(s1.xy_mc(1,:), s1.xy_mc(2,:), 'y');
plot(s1.xy_ekf(1,:), s1.xy_ekf(2,:), 'm');
plot(s1.xy_ukf(1,:), s1.xy_ukf(2,:), 'g');
plot(s1.xy_ckf(1,:), s1.xy_ckf(2,:), 'r');

title('Estimation method comparison');
xlabel('x');
ylabel('y');

test = {'Samples', '$\hat{\mu}_{mc}$', '$\hat{\mu}_{ekf}$', '$\hat{\mu}_{ukf}$', ...
    '$\hat{\mu}_{ckf}$', '$3\sigma_{mc}$', '$3\sigma_{ekf}$', '$3\sigma_{ukf}$', ...
     '$3\sigma_{ckf}$'
    };
legend(test{:}, 'Interpreter','latex', 'FontSize', 12)


% Plot scenario 2
figure(2); clf; hold on; grid on;
%Samples
plot(s2.y(1,:), s2.y(2,:), '.')
%Mean
plot(s2.y_mc(1), s2.y_mc(2), 'yo', options{:});
plot(s2.x_ekf(1), s2.x_ekf(2), 'mo', options{:});
plot(s2.x_ukf(1), s2.x_ukf(2), 'go', options{:});
plot(s2.x_ckf(1), s2.x_ckf(2), 'ro', options{:});
% 3-Sigma Ellipse
plot(s2.xy_mc(1,:), s2.xy_mc(2,:), 'y');
plot(s2.xy_ekf(1,:), s2.xy_ekf(2,:), 'm');
plot(s2.xy_ukf(1,:), s2.xy_ukf(2,:), 'g');
plot(s2.xy_ckf(1,:), s2.xy_ckf(2,:), 'r');

title('Estimation method comparison');
xlabel('x');
ylabel('y');

test = {'Samples', '$\hat{\mu}_{mc}$', '$\hat{\mu}_{ekf}$', '$\hat{\mu}_{ukf}$', ...
    '$\hat{\mu}_{ckf}$', '$3\sigma_{mc}$', '$3\sigma_{ekf}$', '$3\sigma_{ukf}$', ...
     '$3\sigma_{ckf}$'
    };
legend(test{:}, 'Interpreter','latex', 'FontSize', 12)

%% 2 Non-linear Kalman filtering
%% a 
rng(970926)
x0 = [0 0 14 0 0]';
P0 = diag([10 10 2 (pi/180) (5*pi/180)].^2);

s1 = [-200; 100];
s2 = [-200, -100];
Ts = 1;

N = 100;
f = @(x) coordinatedTurnMotion(x, Ts);
h = @(x) dualBearingMeasurement(x, s1, s2);

% Case 1
Q = diag([0 0 1 0 pi/180].^2);
R = diag([10*pi/180 0.5*pi/180].^2);
R = zeros(2);

x = genNonLinearStateSequence(x0, P0, f, Q, N);
y.samples = genNonLinearMeasurementSequence(x, h, R);


figure(1); clf; hold on; grid on;
%Sensors
plot(s1(1), s1(2), 'ro');
plot(s2(1), s2(2), 'ro');
% Measuremsnts
% Calculate intersection
t1 = y.samples(1,:);
t2 = y.samples(2,:);

x_coordinate = (tan(t1)*s1(1) - tan(t2)*s2(1) + s2(2) - s1(2)) ./ (tan(t1) - tan(t2));
y_coordinate = tan(t1).*(x_coordinate - s1(1)) + s1(2);


n = 1;
for type = ["EKF", "UKF", "CKF"]
    figure(n); clf; hold on; grid on;
    
    % True state
    plot(x(1,:), x(2,:));
    % Measurements
    plot(x_coordinate, y_coordinate, 'r.');
    % Filter
    [xf, Pf] = nonLinearKalmanFilter(y.samples, x0, P0, f, Q, h, R, type);
    plot(xf(1,:), xf(2,:));
    
    %for i = 1:5:N
    %    xy = sigmaEllipse2D( xf(1:2,i), Pf(1:2,1:2,i), 3, 500);
    %    plot(xy(1,:), xy(2,:));
    %end
    
    n = n + 1;
end



% True state
plot(x(1,:), x(2,:));
% Filters
%plot(y.mu.ckf(:,i)



legend('Sensor 1', 'Sensor 2', 'Samples', 'True state');

















