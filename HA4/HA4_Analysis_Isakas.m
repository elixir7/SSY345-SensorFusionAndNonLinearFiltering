%% HA4 Analysis Isakas
%% 2 - Analysis
clear; clc; close all;
%% 2.1 Smoothing
clc

rng(6543);
N = 600;
type = "CKF";

%Sensor positions
S = repmat([280; -140],1,N);
   
%Sampling Time
Ts = 0.1;

% Measurement noise is assumed to be known
R = diag([15 2*pi/180].^2);

% Tuned Process noise from HA3
k1 = 0; k2 = 1/4;
Q = diag([0 0 k1*1 0 k2*pi/180].^2);


% Generate true track (state sequence)
omega = zeros(1,N+1);
% Turn rate
omega(200:400) = -pi/201/Ts;
% Initial state
x0 = [0 0 20 0 omega(1)]';
% Allocate memory
x = zeros(length(x0),N+1);
x(:,1) = x0;
% Create true track
for i=2:N+1
    % Simulate
    x(:,i) = coordinatedTurnMotion(x(:,i-1), Ts);
    % Set turn-rate
    x(5,i) = omega(i);
end

% Generate measurements
h = @(x) rangeBearingMeasurements(x, S(:,1));
y = genNonLinearMeasurementSequence(x, h, R);
[x_coord, y_coord] = samplesToState(y, S);


% Initial prior & Covariance for Filter
x0 = [0 0 0 0 0]';
%x0 = [0 0 20 0 0]'; %True initial
P0 = diag([10 10 10 5*pi/180 pi/180].^2);

%Filter with normal Kalman
f = @coordinatedTurnMotion;
h = @rangeBearingMeasurements;
[xs, Ps, xf, Pf, xp, Pp] = nonLinRTSsmoother(y, x0, P0, f, Ts, Q, S, h, R, @sigmaPoints, type);

figure(1); clf; title("Filter vs Smoother");
for i = 1:2
    subplot(1,2,i); hold on; grid on;
    % True track
    plot(x(1,:), x(2,:), '-b', 'LineWidth', 2)
    % Measurements
    plot(x_coord, y_coord, '.m')
    
    if i == 1
        estimator =  type + " Filter";
        plot([x0(1) xf(1,:)], [x0(2) xf(2,:)], '-r', 'LineWidth', 2)
        % 3 Sigma curves
        for k = 1:5:N
            xy = sigmaEllipse2D( xf(1:2,k), Pf(1:2,1:2,k), 3, 500);
            plot(xy(1,:), xy(2,:),  'Color', '#4DBEEE');
        end
    else
        estimator =  "RTS Smoother";
        plot(xs(1,:), xs(2,:), '-r', 'LineWidth', 2)
        % 3 Sigma curves
        for k = 1:5:N
            xy = sigmaEllipse2D(xs(1:2,k), Ps(1:2,1:2,k), 3, 500);
            plot(xy(1,:), xy(2,:),  'Color', '#4DBEEE');
        end
    end

    
    
    %Initial state
    %plot(x0(1), x0(2), 'ro', 'MarkerSize', 15, 'LineWidth', 2);
    %text(x0(1)-50, x0(2)-10, "X0", "FontWeight", "bold")
    
    title(estimator)
    xlabel("X Pos");
    ylabel("Y Pos");
    legend('True state', "Measurements", estimator, '$3\sigma$ curves', 'Interpreter', 'latex')
end





%% b
j = 150;
y(1,j) = y(1, j) + 100;
[x_coord_new, y_coord_new] = samplesToState(y, S);
[xs, Ps, xf, Pf, xp, Pp] = nonLinRTSsmoother(y, x0, P0, f, Ts, Q, S, h, R, @sigmaPoints, type);

figure(2); clf; title("Wrong measurement at k = 150");
for i = 1:2
    subplot(1,2,i); hold on; grid on;
    % True track
    plot(x(1,:), x(2,:), '-b', 'LineWidth', 2)
    % Measurements
    plot(x_coord_new, y_coord_new, '.m')
    
    % Filter and Smoother
    if i == 1
        estimator =  type + " Filter";
        plot([x0(1) xf(1,:)], [x0(2) xf(2,:)], '-r', 'LineWidth', 2)
        % 3 Sigma curves
        for k = 1:5:N
            xy = sigmaEllipse2D( xf(1:2,k), Pf(1:2,1:2,k), 3, 500);
            plot(xy(1,:), xy(2,:),  'Color', '#4DBEEE');
        end
    else
        estimator =  "RTS Smoother";
        plot(xs(1,:), xs(2,:), '-r', 'LineWidth', 2)
        % 3 Sigma curves
        for k = 1:5:N
            xy = sigmaEllipse2D(xs(1:2,k), Ps(1:2,1:2,k), 3, 500);
            plot(xy(1,:), xy(2,:),  'Color', '#4DBEEE');
        end
    end
    
    % Wrong measurement
    plot(x_coord_new(j), y_coord_new(j), 'ro', 'LineWidth', 2);
    text(x_coord_new(j)-200, y_coord_new(j)-10, "Wrong measurement", "FontWeight", "bold")
    
    % Actual measurement
    plot(x_coord(j), y_coord(j),  'ro', 'LineWidth', 2);
    text(x_coord(j), y_coord(j)+15, "Actual measurement", "FontWeight", "bold")
    
    title(estimator)
    xlabel("X Pos");
    ylabel("Y Pos");
    legend('True state', "Measurements", estimator, '$3\sigma$ curves', 'Interpreter', 'latex')
end



%% 2.2 a
clc
rng(970427);

Q = 1.5;
R = 2.5;
x0 = 2;
P0 = 6;

N = 100; %Number of particles

A = 1;
f = @(x) A*x; % Process Model
C = 1; 
h = @(x) C*x; % Measurement Model
T = 0.1; % Sampling time
K = 2/T; % 2s long data

X = genLinearStateSequence(x0, P0, A, Q, K);
Y = genLinearMeasurementSequence(X, C, R);

[X_kal, P_kal] = kalmanFilter(Y, x0, P0, A, Q, C, R);
resampling = false;
kernal = 1;
plotFunc_handle  = @(k, Xk, Xkmin1, Wk, j) ...
                    (plotPostPdf(k, Xk, Wk, X_kal, P_kal, resampling, kernal));
[xfp, Pfp, Xp, Wp] = pfFilter(x0, P0, Y, f, Q, h, R, N, resampling, plotFunc_handle); %Plot func not in arguments
resampling = true;
plotFunc_handle  = @(k, Xk, Xkmin1, Wk, j) ...
                    plotPostPdf(k, Xk, Wk, X_kal, P_kal, resampling, kernal);
[xfp_resamp, Pfp_resamp, Xp_resamp, Wp_resamp] = pfFilter(x0, P0, Y, f, Q, h, R, N, resampling, plotFunc_handle); %Plot func not in arguments



figure(1); clf; hold on; grid on;
plot((0:K).*T ,X, 'b', 'LineWidth', 2)              % True state
plot((1:K).*T,Y, '*m')                              % Measurements
plot((0:K).*T, [x0 X_kal], '-r', 'LineWidth', 2)    % Kalman filter estimate
var = P_kal(:,:,:);                                 % Kalman 3 sigma    
plot((0:K).*T, [x0 X_kal] + 3*sqrt([P0 var(:)']),'--', 'color', '#f07382');
plot((0:K).*T, [x0 X_kal] - 3*sqrt([P0 var(:)']), '--', 'color', '#f07382','HandleVisibility','off');

plot((0:K).*T, [x0 xfp], '-g', 'LineWidth', 2)     % Particle filter estimate without resampling   
var = Pfp(:,:,:);                                   % Particle 3 sigma without resampling
plot((0:K).*T, [x0 xfp] + 3*sqrt([P0 var(:)']),'--', 'color', '#74c971');
plot((0:K).*T, [x0 xfp] - 3*sqrt([P0 var(:)']), '--', 'color', '#74c971','HandleVisibility','off');


plot((0:K).*T, [x0 xfp_resamp], '-m', 'LineWidth', 2)     % Particle filter estimate with resampling  
var = Pfp_resamp(:,:,:);                                   % Parficle 3 sigma with resampling
plot((0:K).*T, [x0 xfp_resamp] + 3*sqrt([P0 var(:)']),'--', 'color', '#d687d2');
plot((0:K).*T, [x0 xfp_resamp] - 3*sqrt([P0 var(:)']), '--', 'color', '#d687d2', 'HandleVisibility','off');

title('Kalman vs Particle filter')
xlabel('Time [s]');
ylabel('State value');
legend('True state ', 'Measurements', 'Kalman','$\pm 3\sigma$ Kalman', ...
    'PF w/o resampling', '$\pm 3\sigma$ PF w/o resampling',  ...
    'PF w/ resampling', '$\pm 3\sigma$ PF w/ resampling', 'Interpreter','latex')



mse.kalman = sum((X(2:end)-X_kal).^2);
mse.pf = sum((X(2:end)-xfp).^2);

%% b/c

figure(1); clf;
subplot(2,1,1); hold on; grid on;
N = 50; %Number of particles
resampling = false;
rng(1997);
[xfp, Pfp, Xp, Wp] = pfFilter(x0, P0, Y, f, Q, h, R, N, resampling, @plotPartTrajs); %Plot func not in arguments

plot(0:K, X, '-r', 'LineWidth', 1)    % True state
xlabel("Sample [k]")
ylabel ("State value")

subplot(2,1,2); hold on; grid on;
N = 50; %Number of particles
resampling = true;
rng(1997);
[xfp, Pfp, Xp, Wp] = pfFilter(x0, P0, Y, f, Q, h, R, N, resampling, @plotPartTrajs); %Plot func not in arguments
plot(0:K, X, '-r', 'LineWidth', 1)    % True state

xlabel("Sample [k]")
ylabel ("State value")





























