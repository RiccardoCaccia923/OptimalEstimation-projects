%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MONTECARLO ANALYSIS FOR LINEAR CALIBRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RMSE vs NOISE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_datasheet = 1.02;
b_datasheet = 0.15;
K_true = [1.02,0,0;0,0.99,0;0,0,1.05];
b_true = [0.13 ; 0.16 ; 0.15];
sigma_vec = linspace(0.001,5,1000);
N_sample = 10; 

g0 = 9.81;
g0_vect = [0;0;g0];

n_meas = 20;
phi = -180 + 360 * randn(1,n_meas); 
theta = -90 + 180 * randn(1,n_meas); 
psi = zeros(1,n_meas);

rmse_final = zeros(length(sigma_vec),1);
rmse_init = zeros(length(sigma_vec),1);

for j = 1:length(sigma_vec)
    sigma = sigma_vec(j);
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
    
    % calibration
    A_calib = inv(K_est_tot) * (a_meas_exp - b_est_tot);
    
    norm_init = vecnorm(a_meas_exp,2,1);
    norm_final = vecnorm(A_calib,2,1);
    
    rmse_init(j) = sqrt(mean((norm_init - g0).^2));
    rmse_final(j) = sqrt(mean((norm_final - g0).^2));

end
    
figure
plot(sigma_vec,rmse_final,'DisplayName','RMSE_{post}')
hold on
plot(sigma_vec,rmse_init,'DisplayName','RMSE_{pre}')
% plot(sigma_vec,sigma_vec,'k','LineWidth',2)
grid on
title('RMSE vs increasing noise')
xlabel('$\sigma$','Interpreter','latex')
ylabel('RMSE','Interpreter','latex')
legend('Location','southeast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RMSE vs # of measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

k_datasheet = 1.02;
b_datasheet = 0.15;
K_true = [1.02,0,0;0,0.99,0;0,0,1.05];
b_true = [0.13 ; 0.16 ; 0.15];
sigma = 0.27;
N_sample = 10; 

g0 = 9.81;
g0_vect = [0;0;g0];

n_vec = (1:1:100);

for j = 1:length(n_vec)
    
    n_meas = n_vec(j);
    phi = -180 + 360 * randn(1,n_meas);  
    theta = -90 + 180 * randn(1,n_meas);
    psi = zeros(1,n_meas);

    A_real = zeros(3,n_meas);
    
    for i = 1:length(theta)
    
        Rx = [1,0,0;0,cos(phi(i)),sin(phi(i));0,-sin(phi(i)),cos(phi(i))];
        Ry = [cos(theta(i)),0,-sin(theta(i));0,1,0;sin(theta(i)),0,cos(theta(i))];
        Rz = [cos(psi(i)),sin(psi(i)),0;-sin(psi(i)),cos(psi(i)),0;0,0,1];
        R = Rx * Ry * Rz;
    
        A_real(:,i) = R * g0_vect;
      
    end
    
    a_meas = K_true * A_real + b_true + sigma * randn(3,n_meas);
    A_real_exp = repelem(A_real,1,N_sample);
    
    noise = sigma * randn(3,n_meas*N_sample);
    a_meas_exp = K_true * A_real_exp + b_true + noise;

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
    
    % calibration
    A_calib = inv(K_est_tot) * (a_meas_exp - b_est_tot);
    
    norm_init = vecnorm(a_meas_exp,2,1);
    norm_final = vecnorm(A_calib,2,1);
    
    rmse_init(j) = sqrt(mean((norm_init - g0).^2));
    rmse_final(j) = sqrt(mean((norm_final - g0).^2));

end
    
figure
plot(n_vec,rmse_final,'LineWidth',1.5,'DisplayName','RMSE_{post}')
hold on
plot(n_vec,rmse_init,'LineWidth',1.5,'DisplayName','RMSE_{pre}')
grid on
title('RMSE vs increasing number of measurements')
xlabel('$\#$ of measurements','Interpreter','latex')
ylabel('RMSE','Interpreter','latex')
legend('Location','southeast')