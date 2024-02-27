function [thetahat, xhat] = lms(x, y, N, mu)
% [thetahat, xhat] = lms(x, y, N, muu)
%
% x         - Observed Data sequence
% y         - Signal sequence
% N         - filter order
% muu       - Step size (learning rate)
% thetahat  - Matrix with estimates of theta. 
% xhat      - Estimate of x
%
% lms: The Least-Mean Square Algorithm
%
% Estimator: xhat(n) = Y^{T}(n)thetahat(n-1)
%            thetahat(n) = thetahat(n − 1) + μY(n){x(n) − xhat(n)}


M = length(x); % Length of data sequences
% Initialize xhat and thetahat
xhat = zeros(M, 1); 
thetahat = zeros(N:M+1); 
% Loop
for n = N:M   
    % Y(n)=[y(n), y(n-1), ...,y(n-N+1)]
    Y_n = y(n:-1:n-N+1);
    % xhat(n) = Y^{T}(n)thetahat(n-1)
    xhat(n) = Y_n' * thetahat(:n-1);  
    % thetahat(n) = thetahat(n − 1) + μY(n){x(n) − xhat(n)} 
    thetahat(:n) = thetahat(:n-1) + mu * Y_n * (x(n) - xhat(n));
end

end
