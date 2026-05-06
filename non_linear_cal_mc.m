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
sigma_vec = linspace(0.001,3,800);
N_sample = 100; 

g0 = 9.81;
g0_vect = [0;0;g0];

n_meas = 20;
phi = -180 + 360 * randn(1,n_meas);
theta = -90 + 180 * randn(1,n_meas); 
psi = zeros(1,n_meas);

rmse_pre_gn = zeros(length(sigma_vec),1);
rmse_post_gn = zeros(length(sigma_vec),1);
iter = zeros(length(sigma_vec),1);

for l = 1:length(sigma_vec)
    
    sigma = sigma_vec(l);

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
    
    a = a_meas_exp;
    N = N_sample*n_meas;
    
    x = [1,1,1,0,0,0]';
    J = zeros(N,6);
    E = zeros(N,1);
    
    max_it = 50;
    toll = 1e-6;
    
    for q = 1 : max_it 
    
        for j = 1:size(J,1)
    
            ax = a(1,j);
            ay = a(2,j);
            az = a(3,j);
    
            Ax = x(1) * ax + x(4);
            Ay = x(2) * ay + x(5);
            Az = x(3) * az + x(6);
        
            A_norm = sqrt(Ax^2+Ay^2+Az^2);
            
            E(j) = A_norm - g0;
        
            J(j,1) = (Ax * ax) / A_norm;
            J(j,2) = (Ay * ay) / A_norm;
            J(j,3) = (Az * az) / A_norm;
            J(j,4) = Ax / A_norm;
            J(j,5) = Ay / A_norm;
            J(j,6) = Az / A_norm;
        
        end
    
        delta_x = (J' * J) \ (J' * E);
        
        x = x - delta_x;
        
        iter(l) = q;
        
        if norm(delta_x) < toll
            % fprintf('Gauss-Newton has converged in %d iterations! \n',iter)
            break;
        elseif norm(delta_x) > toll && q == max_it
            % fprintf('Gauss-Newton cannot converge in %d iterations! \n',max_it)
        end
    
    end
    
    S_est = [x(1),0,0;0,x(2),0;0,0,x(3)];
    o_est = [x(4),x(5),x(6)];
    
    K_eq = inv(S_est);
    b_eq = -1 * o_est / S_est;
    
    A_calib = zeros(3,N);
    for w = 1:N
        A_calib(1,w) = S_est(1,1) * a(1,w) + o_est(1);
        A_calib(2,w) = S_est(2,2) * a(2,w) + o_est(2);
        A_calib(3,w) = S_est(3,3) * a(3,w) + o_est(3);
    end
    
    norm_pre_gn = vecnorm(a,2,1);
    norm_post_gn = vecnorm(A_calib,2,1);
    
    rmse_pre_gn(l) = sqrt(mean((norm_pre_gn-g0).^2));
    rmse_post_gn(l) = sqrt(mean((norm_post_gn-g0).^2));

end

figure
plot(sigma_vec,rmse_post_gn,'DisplayName','RMSE_{postGN}')
hold on
plot(sigma_vec,rmse_pre_gn,'DisplayName','RMSE_{preGN}')
grid on
%title('Gauss-Newton RMSE vs increasing Noise')
xlabel('$\sigma$','Interpreter','latex')
ylabel('Residual Errors $m/s^2$','Interpreter','latex')
legend show

figure
plot(sigma_vec,iter,'DisplayName','GN_{iter}')
grid on
%title('Gauss-Newton iterations vs increasing Noise')
xlabel('$\sigma$','Interpreter','latex')
ylabel('$\#$ of iterations','Interpreter','latex')
legend show

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

rmse_pre_gn = zeros(length(n_vec),1);
rmse_post_gn = zeros(length(n_vec),1);
iter = zeros(length(n_vec),1);

for l = 1:length(n_vec)
    
    n_meas = n_vec(l);

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
    
    a = a_meas_exp;
    N = N_sample*n_meas;
    
    x = [1,1,1,0,0,0]';
    J = zeros(N,6);
    E = zeros(N,1);
    
    max_it = 50;
    toll = 1e-6;
    
    for q = 1 : max_it 
    
        for j = 1:size(J,1)
    
            ax = a(1,j);
            ay = a(2,j);
            az = a(3,j);
    
            Ax = x(1) * ax + x(4);
            Ay = x(2) * ay + x(5);
            Az = x(3) * az + x(6);
        
            A_norm = sqrt(Ax^2+Ay^2+Az^2);
            
            E(j) = A_norm - g0;
        
            J(j,1) = (Ax * ax) / A_norm;
            J(j,2) = (Ay * ay) / A_norm;
            J(j,3) = (Az * az) / A_norm;
            J(j,4) = Ax / A_norm;
            J(j,5) = Ay / A_norm;
            J(j,6) = Az / A_norm;
        
        end
    
        delta_x = (J' * J) \ (J' * E);
        
        x = x - delta_x;
        
        iter(l) = q;
        
        if norm(delta_x) < toll
            % fprintf('Gauss-Newton has converged in %d iterations! \n',iter)
            break;
        elseif norm(delta_x) > toll && q == max_it
            % fprintf('Gauss-Newton cannot converge in %d iterations! \n',max_it)
        end
    
    end
    
    S_est = [x(1),0,0;0,x(2),0;0,0,x(3)];
    o_est = [x(4),x(5),x(6)];
    
    K_eq = inv(S_est);
    b_eq = -1 * o_est / S_est;
    
    A_calib = zeros(3,N);
    for w = 1:N
        A_calib(1,w) = S_est(1,1) * a(1,w) + o_est(1);
        A_calib(2,w) = S_est(2,2) * a(2,w) + o_est(2);
        A_calib(3,w) = S_est(3,3) * a(3,w) + o_est(3);
    end
    
    norm_pre_gn = vecnorm(a,2,1);
    norm_post_gn = vecnorm(A_calib,2,1);
    
    rmse_pre_gn(l) = sqrt(mean((norm_pre_gn-g0).^2));
    rmse_post_gn(l) = sqrt(mean((norm_post_gn-g0).^2));

end

figure
plot(1:n_meas,rmse_post_gn,'LineWidth',1.5,'DisplayName','RMSE_{postGN}')
hold on
plot(1:n_meas,rmse_pre_gn,'LineWidth',1.5,'DisplayName','RMSE_{preGN}')
grid on
%title('Gauss-Newton RMSE vs $\#$ of measurements','Interpreter','latex')
xlabel('$\#$ of measurements','Interpreter','latex')
ylabel('Residual Errors $m/s^2$','Interpreter','latex')
legend show

figure
plot(1:n_meas,iter,'LineWidth',1.5,'DisplayName','GN_{iter}')
grid on
%title('Gauss-Newton iterations vs $\#$ of measurements','Interpreter','latex')
xlabel('$\#$ of measurements','Interpreter','latex')
ylabel('$\#$ of iterations','Interpreter','latex')
legend show


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RMSE vs Initial Guess
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

k_datasheet = 1.02;
b_datasheet = 0.15;
K_true = [1.02,0,0;0,0.99,0;0,0,1.05];
b_true = [0.13 ; 0.16 ; 0.15];
sigma = 0.027;
N_sample = 100; 

g0 = 9.81;
g0_vect = [0;0;g0];

n_meas = 20;
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



a = a_meas_exp;
N = N_sample*n_meas;

J = zeros(N,6);
E = zeros(N,1);

max_it = 50;
toll = 1e-6;

n_sim = 100;

rmse_pre_gn = zeros(n_sim,1);
rmse_post_gn = zeros(n_sim,1);
iter = zeros(n_sim,1);

for l = 1 : n_sim

    max_err_s = 0.3;
    max_err_o = 0.3;
    
    rand_s = (rand(1,3) - 0.5) * 2;
    rand_o = (rand(1,3) - 0.5) * 2;

    s_var = 1 + max_err_s * rand_s;
    o_var = 0 + max_err_o * rand_o;
    
    x = [s_var,o_var]';
    
    for q = 1 : max_it 
    
        for j = 1:size(J,1)
    
            ax = a(1,j);
            ay = a(2,j);
            az = a(3,j);
    
            Ax = x(1) * ax + x(4);
            Ay = x(2) * ay + x(5);
            Az = x(3) * az + x(6);
        
            A_norm = sqrt(Ax^2+Ay^2+Az^2);
            
            E(j) = A_norm - g0;
        
            J(j,1) = (Ax * ax) / A_norm;
            J(j,2) = (Ay * ay) / A_norm;
            J(j,3) = (Az * az) / A_norm;
            J(j,4) = Ax / A_norm;
            J(j,5) = Ay / A_norm;
            J(j,6) = Az / A_norm;
        
        end
    
        delta_x = (J' * J) \ (J' * E);
        
        x = x - delta_x;
        
        iter(l) = q;
        
        if norm(delta_x) < toll
            %fprintf('Gauss-Newton has converged in %d iterations! \n',iter)
            break;
        elseif norm(delta_x) > toll && q == max_it
            %fprintf('Gauss-Newton cannot converge in %d iterations! \n',max_it)
        end
    
    end
    
    S_est = [x(1),0,0;0,x(2),0;0,0,x(3)];
    o_est = [x(4),x(5),x(6)];
    
    K_eq = inv(S_est);
    b_eq = -1 * o_est / S_est;

    A_calib = zeros(3,N);

    for w = 1:N
        A_calib(1,w) = S_est(1,1) * a(1,w) + o_est(1);
        A_calib(2,w) = S_est(2,2) * a(2,w) + o_est(2);
        A_calib(3,w) = S_est(3,3) * a(3,w) + o_est(3);
    end
    
    norm_pre_gn = vecnorm(a,2,1);
    norm_post_gn = vecnorm(A_calib,2,1);
    
    rmse_pre_gn(l) = sqrt(mean((norm_pre_gn-g0).^2));
    rmse_post_gn(l) = sqrt(mean((norm_post_gn-g0).^2));

end

figure
plot(1:n_sim,rmse_post_gn,'LineWidth',1.5,'DisplayName','RMSE_{postGN}')
hold on
%plot(1:n_sim,rmse_pre_gn,'LineWidth',1.5,'DisplayName','RMSE_{preGN}')
grid on
ylim=[-1,1];
%title('Gauss-Newton RMSE vs $\#$ of measurements','Interpreter','latex')
xlabel('time','Interpreter','latex')
ylabel('Residual Errors $m/s^2$','Interpreter','latex')
legend show

figure
plot(1:n_sim,iter,'LineWidth',1.5,'DisplayName','GN_{iter}')
grid on
%title('Gauss-Newton iterations vs $\#$ of measurements','Interpreter','latex')
xlabel('time','Interpreter','latex')
ylabel('$\#$ of iterations','Interpreter','latex')
legend show