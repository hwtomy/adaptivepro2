% Project2 of Tiou Wang and Liping Bai

% Problem Setup and Notation System
% 

% Optimal Linear Estimator 
% ΣYY(n)θopt(n)=ΣYx(n) 
% has a computational cost of O(n^3)

% The Steepest-Descent Algorithm derived with Expected MSE
% MSE(n,theta) = rxx(n,n) − 2theta^{T}ΣYx(n) + theta^{T}ΣYY(n)theta
% theta_hat(n) = theta_hat(n−1) − μ/2*∂MSE(n,theta)/∂theta∣theta=theta_hat(n−1)
% theta_hat(n) = theta_hat(n − 1) − μ{ΣYY(n)theta_hat(n − 1) − ΣYx(n)}
% has a computational cost of O(n^2)

% The Least Mean Square (LMS) Algorithm derived with instantaneous error
% MSE_hat(n, theta) = {x(n) − theta^{T}Y(n)}^2
% theta_hat(n) = theta_hat(n−1)−μ/2*∂MSE_hat(n,theta)/∂theta∣theta=theta_hat(n−1)
% theta_hat(n) = theta_hat(n − 1) + μY(n)*{x(n) − x_hat(n)}
% x_hat(n)=Y^{T}(n)theta_hat(n-1)
% has a computational cost of O(n)

% The Normalized LMS (NLMS) Algorithm has normalized step size
% μ(n) = μ_bar/{c + ‖Y (n)‖_2Norm}

% The Recursive Least Square (RLS) Algorithm
% Σ_hat_RLS_YY(n) = Y(n)Y^T(n) + λΣ_hat_RLS_YY(n − 1)
% Σ_hat_RLS_Yx(n) = Y(n)x(n) + λΣ_hat_RLS_Yx(n − 1)
% θ_hat_(n) = ˆθ(n − 1) +μInv{Σ_hat_RLS_YY(n)}
%                       {Σ_hat_RLS_Yx(n)−Σ_hat_RLS_YY(n)θ_hat(n − 1)}
% has a computational cost of O(n^3)


clc;
clear all;
clc;

[s,fs]=audioread("EQ2401Project2data2024.wav");
[s_additional, fs]=audioread("EQ2401Project2_bonus_task2024.wav");
% plot(periodogram(s,hanning(length(s)), fs));

D = 0;
delay = dsp.Delay(D);
x = delay(s);

%%lms
mu = 0.01; 
order = 256;
ga = 0.1;
w = lms(s,x,mu,order);
sn=transpose(s(1:order-1));
for i=order:length(s)
    sn1 = s(i:-1:i-order+1)'*w(:,i-1);
    sn=[sn,sn1];
end
sn = sn';
sv_lms = s-sn;
% soundsc(sv,fs);

%%rls
mu = 0.6; 
order = 512;
w = rls(s,x,mu,order);
sn=transpose(s(1:order-1));
for i=order:length(s)
    sn1 = s(i:-1:i-order+1)'*w(:,i-1);
    sn=[sn,sn1];
end
sn = sn';
sv_rls = s-sn;
% soundsc(sv,fs);

%%nlms
mu = 0.999; 
order = 512;
w = nlms(s,x,mu,order);
sn=transpose(s(1:order-1));
for i=order:length(s)
    sn1 = s(i:-1:i-order+1)'*w(:,i-1);
    sn=[sn,sn1];
end
sn = sn';
sv_nlms = s-sn;
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
w = adam(s,x,mu,order,b1,b2);
sn=transpose(s(1:order-1));
for i=order:length(s)
    sn1 = s(i:-1:i-order+1)'*w(:,i-1);
    sn=[sn,sn1];
end
sn = sn';
sv_adam = s-sn;
%soundsc(sv,fs);


%% Plotting
figure; % Create a new figure

plot(s, 'LineWidth', 1.5); % Plot original signal 
hold on; % Hold the current plot
plot(sv_lms); 
plot(sv_rls);
plot(sv_nlms);
plot(sv_adam);

xlabel('Sample');
ylabel('Amplitude');
title('Comparison of Adaptive Filter Outputs');
legend('Original Signal', 'LMS Output', 'RLS Output', 'NLMS Output'); 
grid on; 





