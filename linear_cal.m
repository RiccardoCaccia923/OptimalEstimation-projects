%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ACCELEROMETER LINEAR CALIBRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear
clc

%% DATA
% accelerometer datasheet
k_datasheet = 1.02;
b_datasheet = 0.15;
K_true = [1.02,0,0;0,0.99,0;0,0,1.05];
b_true = [0.13 ; 0.16 ; 0.15];
sigma = 0.027;
N_sample = 100; % sampling frequency

% gravity local vector
g0 = 9.81;
g0_vect = [0;0;g0];

% accelerometer measurements
n_meas = 20;
phi = -180 + 360 * randn(1,n_meas);  % random between -180 and 180
theta = -90 + 180 * randn(1,n_meas); % random between -90 and 90
psi = zeros(1,n_meas);

%% TRUE ACCELERATIONS AND MEASUREMENT
A_real = zeros(3,n_meas);

for i = 1:length(theta)

    % define rotation matrix
    Rx = [1,0,0;0,cos(phi(i)),sin(phi(i));0,-sin(phi(i)),cos(phi(i))];
    Ry = [cos(theta(i)),0,-sin(theta(i));0,1,0;sin(theta(i)),0,cos(theta(i))];
    Rz = [cos(psi(i)),sin(psi(i)),0;-sin(psi(i)),cos(psi(i)),0;0,0,1];
    R = Rx * Ry * Rz;

    % compute true acceleration for each measurement
    A_real(:,i) = R * g0_vect;
  
end

a_meas = K_true * A_real + b_true + sigma * randn(3,n_meas);
A_real_exp = repelem(A_real,1,N_sample);

% compute measured accelerations 
noise = sigma * randn(3,n_meas*N_sample);
a_meas_exp = K_true * A_real_exp + b_true + noise;


% acc measurement plot
time = 1:n_meas*N_sample;
figure
plot(time,a_meas_exp(1,:),'r','LineWidth',1.5,'DisplayName','a_x')
hold on
plot(time,a_meas_exp(2,:),'b','LineWidth',1.5,'DisplayName','a_y')
plot(time,a_meas_exp(3,:),'g','LineWidth',1.5,'DisplayName','a_z')
grid on
title('Accelerometer Measurements')
xlabel('Time [s]')
ylabel('Measured Acceleration [m/s^2]')
legend

%% OLS Algorithm
K_est_tot = zeros(3);
b_est_tot = zeros(3,1);
residual_tot = zeros(n_meas*N_sample,3);

for i = 1:3
    Y = a_meas_exp(i,:)';
    H = [A_real_exp(i,:)', ones(n_meas*N_sample,1)];
    x_est = inv( H' * H ) * H' * Y;
    K_est_tot(i,i) = x_est(1);
    b_est_tot(i) = x_est(2);
    residual_tot(:,i) = Y - H * x_est;
end

K_est_mean = zeros(3);
b_est_mean = zeros(3,1);
residual_mean = zeros(n_meas,3);

for i = 1:3
    Y = a_meas(i,:)';
    H = [A_real(i,:)', ones(n_meas,1)];
    x_est = inv( H' * H ) * H' * Y;
    K_est_mean(i,i) = x_est(1);
    b_est_mean(i) = x_est(2);
    residual_mean(:,i) = Y - H * x_est;
end

% rmse computation
A_calib = inv(K_est_tot) * (a_meas_exp - b_est_tot);

norm_init = vecnorm(a_meas_exp,2,1);
norm_final = vecnorm(A_calib,2,1);

rmse_init = sqrt(mean((norm_init - g0).^2));
rmse_final = sqrt(mean((norm_final - g0).^2));

%% RESULTS
% estimated paraeters
disp('--- OLS ALGORITHM RESULTS ---')
disp('- Scale Factors (K) -----------');
disp('True:'); disp(K_true);
disp('Estimated:'); disp(K_est_tot);
disp('- Bias (b) [m/s^2] ------------');
disp('True:'); disp(b_true);
disp('Estimated:'); disp(b_est_tot);
disp('---------------------------------------------')

disp('--- OLS ALGORITHM RESULTS ---')
disp('- Scale Factors (K) -----------');
disp('True:'); disp(K_true);
disp('Estimated:'); disp(K_est_mean);
disp('- Bias (b) [m/s^2] ------------');
disp('True:'); disp(b_true);
disp('Estimated:'); disp(b_est_mean);

% rmse
disp('- RMSE Analysis ---------------');
fprintf('RMSE Pre-Calibration: %.4f m/s^2 \n', rmse_init)
fprintf('RMSE Post-Calibration: %.4f m/s^2 \n', rmse_final);

% residuals
figure
plot(time,residual_tot(:,1),'r','DisplayName','res_x')
hold on
plot(time,residual_tot(:,2),'b','DisplayName','res_y')
plot(time,residual_tot(:,3),'g','DisplayName','res_z')
grid on
title('OLS Total Residuals')
xlabel('time [s]')
ylabel('residuals')
legend show

figure
plot(1:n_meas,residual_mean(:,1),'ro','LineWidth',1.5,'DisplayName','res_x')
hold on
plot(1:n_meas,residual_mean(:,2),'bo','LineWidth',1.5,'DisplayName','res_y')
plot(1:n_meas,residual_mean(:,3),'go','LineWidth',1.5,'DisplayName','res_z')
grid on
title('OLS Mean Residuals')
xlabel('time [s]')
ylabel('residuals')
legend show
disp('---------------------------------------------')
