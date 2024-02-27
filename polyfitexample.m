% fitting a polynomial by RLS from streaming data 
% MJ 190218
%

%% Generate data
clc,clear,close all

N = 100;
sig2 = .05; %noise variance
x = randn(N,1); %random sampling points on the x-axis

X = [ones(N,1) x x.^2]; % regressor matrix, 2nd order poly
m = size(X,2); %order of polynomial
a = rand(m,1); % true poly coeffs

e = sqrt(sig2)*randn(N,1);
y = X*a + e; % noisy polynomial observations

%sort to plot true graph
[xsort,ind] = sort(x);
Xsort = X(ind,:);
figure(1),movegui
plot(xsort,Xsort*a,'b'),hold on,
title(['Polynomial: y= ' num2str(a(1),2) '+' ...
    num2str(a(2),2) 'x +' num2str(a(3),2) 'x^2'  ])

%% plot for poly coeffs
figure(2), movegui
plot(1:N,a*ones(1,N),'b'),hold on
title('Polynomial coefficients')

%% RLS
% initialize RLS
xhatrls = zeros(N,1);
thrls = zeros(N+1,m);
Prls = 1e2*eye(m);
lambda = 1;

for n=1:N

% rls
  xhatrls(n) = thrls(n,:)*X(n,:)';
  K = Prls*X(n,:)'/(lambda+X(n,:)*Prls*X(n,:)');
  Prls = 1/lambda*(Prls-K*X(n,:)*Prls);  
  thrls(n+1,:) = thrls(n,:)+K'*(y(n)-xhatrls(n));

  figure(1)
  %plot(x(n),y(n),'+'),plot(x(n),xhatrls(n),'or'),pause(.1),drawnow
  plot(x(n),y(n),'+'),pause(.1),drawnow
hold on
  figure(2)
  plot(n+1,thrls(n+1,:),'o'),drawnow

end
hold off

%% plot the polynomial estimated by RLS at the final point
figure(1),plot(xsort,Xsort*thrls(end,:)','r'),hold off

%% compare the batch least squares estimate with that of RLS
als = X\y
thrls(end,:)'

%% 
if 1
% X nonstationary!?
%X'*X/N
(Xsort'*Xsort)/N
(Xsort(1:floor(N/2),:)'*Xsort(1:floor(N/2),:))/(N/2)
(Xsort(floor(N/2)+1:N,:)'*Xsort(floor(N/2)+1:N,:))/(N/2)
end
