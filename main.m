% Project2 of Tiou Wang and Liping Bai

% System Model and Notation System
% Measured signal: z(n) = s(n) + x(n)
% where s(n) is the desired signal and x(n) is correlated with measurable disturbance y(n)
% Estimation of desired signal: shat(n) = z(n) − xhat(n)
% in this project we used Adaptive line Enhancer(ALE) to predict the noise
% i.e. the desired signal is set to be a delayed version of the measured data.

% Optimal Linear Minumum Mean Square Error (LMMSE) Estimator of a statinary process
% ΣYYtheta_opt=ΣYX 
% When the underlying system changes with time n
% MSE(n,theta) = E{(x(n)− x_hat(n))^2} where x_hat(n)=Y^{T}(n)theta
% has a computational cost of O(N^3)

% The Steepest-Descent Algorithm derived with Expected MSE
% MSE(n,theta) = rxx(n,n) − 2theta^{T}ΣYx(n) + theta^{T}ΣYY(n)theta
% theta_hat(n) = theta_hat(n−1) − μ/2*∂MSE(n,theta)/∂theta∣theta=theta_hat(n−1)
% theta_hat(n) = theta_hat(n − 1) − μ{ΣYY(n)theta_hat(n − 1) − ΣYx(n)}
% has a computational cost of O(N^2)

% The Least Mean Square (LMS) Algorithm derived with instantaneous error
% MSE_hat(n, theta) = {x(n) − theta^{T}Y(n)}^2
% theta_hat(n) = theta_hat(n−1)−μ/2*∂MSE_hat(n,theta)/∂theta∣theta=theta_hat(n−1)
% theta_hat(n) = theta_hat(n − 1) + μY(n)*{x(n) − x_hat(n)}
% x_hat(n)=Y^{T}(n)theta_hat(n-1)
% has a computational cost of O(N)

% The Normalized LMS (NLMS) Algorithm has normalized step size
% μ(n) = μ_bar/{c + ‖Y (n)‖_2Norm}

% The Recursive Least Square (RLS) Algorithm
% Σ_hat_RLS_YY(n) = Y(n)Y^T(n) + λΣ_hat_RLS_YY(n − 1)
% Σ_hat_RLS_Yx(n) = Y(n)x(n) + λΣ_hat_RLS_Yx(n − 1)
% θ_hat_(n) = ˆθ(n − 1) +μInv{Σ_hat_RLS_YY(n)}
%                       {Σ_hat_RLS_Yx(n)−Σ_hat_RLS_YY(n)θ_hat(n − 1)}
% has a computational cost of O(N^3)


clc;
clear all;
clc;

[z,fs]=audioread("EQ2401Project2data2024.wav");

% set the filter parameters
mu_lms = 0.001;
mu_sd = 0.001; 
N = 300;
M = length(z);
xhat_lms = zeros(M, 1);
xhat_sd = zeros(M,1); 
thetahat_lms = zeros(N,M+1);
thetahat_sd =  zeros(N,M+1);
delay = 0; 

for n = N+delay:M
    % Y (n) = [y(n), y(n − 1), . . . , y(n − N + 1)]T 
    Y_n = z(n:-1:n-N+1);
    % xhat(n) = Y^{T}(n)thetahat(n-1) 
    xhat_lms(n) = Y_n' * thetahat_lms(:, n-1);  
    % thetahat(n) = thetahat(n − 1) + μY(n){x(n) − xhat(n)}
    thetahat_lms(:, n) = thetahat_lms(:, n-1) + mu_lms * Y_n * (z(n-delay) - xhat_lms(n));

    %Steepest descent
    % theta_hat(n) = theta_hat(n − 1) − μ{ΣYY(n)theta_hat(n − 1) − ΣYx(n)}
    SigmaYY = Y_n * Y_n'; % this is a N by N matrix
    SigmaYx = Y_n * z(n-delay); 
    xhat_sd(n) = Y_n'*thetahat_sd(:,n-1);
    thetahat_sd(:,n) = thetahat_sd(:,n-1) - mu_sd*(SigmaYY*thetahat_sd(:,n-1) - SigmaYx); 
end

s_lms = z-xhat_lms;
s_sd = z-xhat_sd;
soundsc(s_sd, fs);
audiowrite('lms.wav', s_lms,fs);


% Find the index of the most significant coefficient for both LMS and SD
[~, index_lms] = max(abs(thetahat_lms(:, end)));
[~, index_sd] = max(abs(thetahat_sd(:, end)));

% Plot time evolution of the most significant coefficient
figure; % Create a new figure for the plot
hold on; % Allow multiple plots on the same figure

% Plot for LMS
plot(thetahat_lms(index_lms, :), 'DisplayName', sprintf('LMS Coeff %d', index_lms), LineWidth=2); 

% Plot for Steepest Descent
plot(thetahat_sd(index_sd, :), 'DisplayName', sprintf('SD Coeff %d', index_sd)); 

hold off; % Stop overplotting
xlabel('Time (Sample)');
ylabel('Theta Value');
title('Comparison of the Most Significant Theta Coefficient for LMS and Steepest Descent');
legend('show');


%%rls
mu = 0.6; 
order = 512;
w = rls(z,x,mu,order);
sn=transpose(z(1:order-1));
for i=order:length(s)
    sn1 = z(i:-1:i-order+1)'*w(:,i-1);
    sn=[sn,sn1];
end
sn = sn';
sv_rls =z-sn;
% soundsc(sv,fs);

%%nlms
mu = 0.999; 
order = 512;
w = nlms(z,x,mu,order);
sn=transpose(z(1:order-1));
for i=order:length(s)
    sn1 = z(i:-1:i-order+1)'*w(:,i-1);
    sn=[sn,sn1];
end
sn = sn';
sv_nlms = z-sn;
% soundsc(sv,fs);

%%adagradelsm
% mu = 0.05; 
% order = 256;
% w = adag(s,x,mu,order);
% sn=transpose(s(1:order-1));
% for i=order:length(s)
%     sn1 = s(i:-1:i-order+1)'*w(:,i-1);
%     sn=[sn,sn1];
% end
% sn = sn';
% sv = s-sn;
% soundsc(sv,fs);

%%adamlsm
mu = 0.8; 
b1 = 0.2;
b2 = 0.0000001;
order = 256;
w = adam(z,x,mu,order,b1,b2);
sn=transpose(z(1:order-1));
for i=order:length(s)
    sn1 = z(i:-1:i-order+1)'*w(:,i-1);
    sn=[sn,sn1];
end
sn = sn';
sv_adam = z-sn;
%soundsc(sv,fs);


%%Plotting
figure; % Create a new figure

plot(z, 'LineWidth', 1.5); % Plot original signal 
hold on; % Hold the current plot
plot(sv_lms); 
plot(sv_rls);
plot(sv_nlms);
%plot(sv_adam);

xlabel('Sample');
ylabel('Amplitude');
title('Comparison of Adaptive Filter Outputs');
legend('Original Signal', 'LMS Output', 'RLS Output', 'NLMS Output', 'Adam Output'); 
grid on; 