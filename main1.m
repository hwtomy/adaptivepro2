%
%  RLS_LMS_comp:
%     Examples for the slides in Lecture 10/11.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Example 1
n=0:3000;
y=sin(2*pi*n/300);
eta0=0.5;
e=randn(size(n));
x=eta0*y+e;

[etahatLMS,xhatLMS]=lms_vect(x,y,0.01);
[etahatRLS,xhatRLS]=rls_vect(x,y,0.995);
[etahatRLS1,xhatRLS1]=rls_vect(x,y,1);
%[etahatLMS,xhatLMS]=lms(x,y,1,0.01);
%[etahatRLS,xhatRLS]=rls(x,y,1,0.995);
%[etahatRLS1,xhatRLS1]=rls(x,y,1,1);

figure(1)
plot(n,etahatLMS, n,etahatRLS, 'r--', n, etahatRLS1, '-.g')
ylim([-.6,.8])
legend('etahat(LMS) \mu=0.01', 'etahat(RLS) \lambda=0.995', ...
       'etahat(RLS) \lambda=1')
text(500,0,'2/(1-0.995)=400')
text(500,-0.2,'1-0.995=0.005')
text(500,-0.4,'\mu \sigma_{yy}^2=0.01*0.5=0.005')
title('Comparison of RLS and LMS')

%% Example 2
M=12000;
n=1:M;
theta=[.5*ones(1,M/3),1*ones(1,M/3),.5*ones(1,M/3);0.5*ones(1,M)];
%theta=[.5*ones(1,M/3),.5*ones(1,M/3),.5*ones(1,M/3);0.5*ones(1,M)];
SigmaYY=[0.5 1; 1 3];
Y=chol(SigmaYY)'*randn(2,M);
x=sum(theta.*Y)+0.1*randn(1,M);

[thetahatLMS,xhatLMS]=lms_vect(x,Y,0.01);
[thetahatRLS,xhatRLS]=rls_vect(x,Y,0.997);
[thetahatRLS1,xhatRLS1]=rls_vect(x,Y,1);

figure(2)
subplot(2,1,1)
h1=plot(n, theta(1,:), 'k', n,thetahatLMS(:,1), ':b', ...
     n,thetahatRLS(:,1), 'r--', n, thetahatRLS1(:,1), '-.g');
ylim([0,1.5])
xlabel('n'),ylabel('\theta_{1}')
title('Comparison of RLS and LMS')

subplot(2,1,2)
h2=plot(n, theta(2,:), 'k', n,thetahatLMS(:,2), ':b', ...
     n,thetahatRLS(:,2), 'r--', n, thetahatRLS1(:,2), '-.g');
 xlabel('n'),ylabel('\theta_{2}')
legend('True value', 'LMS \mu=0.01', 'RLS \lambda=0.997', ...
       'RLS \lambda=1')
set(h1,'linewidth',2)
set(h2,'linewidth',2)