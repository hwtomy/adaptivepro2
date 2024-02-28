% Project2 of Tiou Wang and Liping Bai

% System Model and Notation System
% Measured signal: z(n) = s(n) + x(n)
% Estimation of desired signal: shat(n) = z(n) − xhat(n)
% in this project we used Adaptive line Enhancer(ALE) to predict the noise
% i.e. the desired signal is set to be a delayed version of the measured data.

% This file is standalone file. Can be excuted without any supporting function.

clc;
clear all;
clc;
% Read data
[z,fs]=audioread("EQ2401Project2data2024.wav");
%[z,fs]=audioread("EQ2401Project2_bonus_task2024.wav");

% set the filter parameters
mu_lms = 0.004;
mu_bar = 0.25;
N = 300;
M = length(z);
D = 0;
delay = dsp.Delay(D);
x= delay(z);
P = .05*eye(N);
lambda = 0.9999999;
c = 1e-10;

%Initiate data structure
xhat_lms = zeros(M, 1);
xhat_nlms = zeros(M, 1);
xhat_rls = zeros(M,1);
%from i to N-1, there is no filtering
for i = 1:N-1
    xhat_lms(i)=z(i);
    xhat_nlms(i)=z(i);
    xhat_rls(i)=z(i);
end
mu_nlms = zeros(M,1);
thetahat_lms = zeros(N,M+1);
thetahat_nlms = zeros(N,M+1);
thetahat_rls =  zeros(N,M+1);
mean_y = zeros(1, M);
var_y = zeros(1, M);

% Filtering
for n = N:M
    % Y (n) = [y(n), y(n − 1), . . . , y(n − N + 1)]T 
    Y_n = z(n:-1:n-N+1);
    % compute for the statistical characteristics of Y_n for plotting
    mean_y(n) = mean(Y_n);
    var_y(n) = var(Y_n);
    
    %%  Least Mean Square (LMS) Algorithm
    % Instantaneous Error MSE(n, theta) = {x(n) − theta^{T}Y(n)}^2
    % theta_hat(n) = theta_hat(n−1)−μ/2*∂MSE_hat(n,theta)/∂theta∣theta=theta_hat(n−1)
    % theta_hat(n) = theta_hat(n − 1) + μY(n)*{x(n) − x_hat(n)}
    % x_hat(n)=Y^{T}(n)theta_hat(n-1)
    xhat_lms(n) = Y_n' * thetahat_lms(:, n-1);  
    thetahat_lms(:, n) = thetahat_lms(:, n-1) + mu_lms * Y_n * (x(n) - xhat_lms(n));

    %% Normalized LMS (NLMS) Algorithm
    xhat_nlms(n) = Y_n' * thetahat_lms(:, n-1); 
    % μ(n) = μ_bar/{c + ‖Y(n)‖L2Norm}
    mu_nlms(n) = mu_bar/(c+Y_n'*Y_n);
    thetahat_nlms(:,n) = thetahat_nlms(:,n-1) + mu_nlms(n) *Y_n*(x(n)- xhat_nlms(n)); 

    %% Recursive Least Square (RLS) Algorithm
   
    % Σ_hat_YY(n) = Y(n)Y^T(n) + λΣ_hat_YY(n − 1)
    % Σ_hat_Yx(n) = Y(n)x(n) + λΣ_hat_Yx(n − 1)
    % θ_hat(n) = ˆθ(n − 1) +μInv{Σ_hat_YY(n)}
    %                       {Σ_hat_Yx(n)−Σ_hat_YY(n)θ_hat(n − 1)}
  
    xhat_rls(n) = Y_n' * thetahat_rls(:,n-1);
    K=P*Y_n/(lambda+Y_n'*P*Y_n);
    P=1/lambda*(P-K*Y_n'*P);
    thetahat_rls(:,n)=thetahat_rls(:,n-1)+K*(x(n)-xhat_rls(n));
end

% Store data
shat_lms = z-xhat_lms;
shat_nlms = z-xhat_nlms;
shat_rls = z-xhat_rls;
soundsc(shat_nlms, fs);
audiowrite('lms.wav', shat_lms,fs);
audiowrite('nlms.wav', shat_nlms,fs);
audiowrite('rls.wav', shat_rls,fs);

% Plot
figure; 
hold on; 
for i=1:10
    % Plot for LMS
    plot(thetahat_lms(i, :), 'Color','red'); 
    % Plot for NLMS
    plot(thetahat_nlms(i, :), 'Color','green'); 
    % Plot for RLS
    plot(thetahat_rls(i, :), 'Color','blue'); 
end
hLMS = plot(NaN,NaN,'Color','red', 'LineWidth',2);
hNLMS = plot(NaN,NaN,'Color','green');
hRLS = plot(NaN,NaN,'Color','blue');
legend([hLMS, hNLMS, hRLS], {'LMS', 'NLMS', 'RLS'});
hold off;
xlabel('Time (Sample)');
ylabel('Theta Value');
title('Comparison of the Most Significant Theta Coefficient for LMS, NLMS, and RLS');

figure; 
%plot(z, 'LineWidth', 1.5);
hold on;
plot(shat_lms, 'LineWidth', 3);
plot(shat_nlms);
plot(shat_rls);
xlabel('Time (Sample)');
ylabel('Amplitude');
title('Comparison of Adaptive Filter Outputs');
legend('LMS Output', 'NLMS Output', 'RLS Output'); 
grid on; 

figure;
plot(mu_nlms);
xlabel('Time (Sample)');
ylabel('μ(n) = μ_bar/{c + ‖Y(n)‖L2Norm}')
grid on;


figure;
subplot(2,1,1); % Plot mean of Y_n
plot(N:M, mean_y(N:M));
title('Mean of y_n Over Time');
xlabel('Time (Sample)');
ylabel('Mean');

subplot(2,1,2); % Plot variance of Y_n
plot(N:M, var_y(N:M));
title('Variance of y_n Over Time');
xlabel('Time (Sample)');
ylabel('Variance');