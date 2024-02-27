function [thetahat,xhat]=lms_vect(x,Y,muu)

    % [thetahat,xhat]=lms_vect(x,Y,muu)
    %
    %	x			- Data sequence
    %	Y			- Data sequence
    %	muu			- Step size
    %	thetahat		- Matrix with estimates of theta. 
    %				  Row n corresponds to the estimate thetahat(n)'
    %	xhat			- Estimate of x
    %
    %
    %
    %  lms: The Least-Mean Square Algorithm
    %
    % 	Estimator: xhat(n)=Y(:,n)^{T} thetahat(n-1)
    %
    %	thetahat is estimated using LMS. 
    %
    %     
    %     Author: 
    %
    
    % Initialize xhat and thetahat
    [N,M]=size(Y);
    xhat=zeros(M,1);
    thetahat=zeros(M,N);
    
    % Loop
    for n=1:M,
      % Generate Y. Set elements of Y that does not exist to zero
      % Estimate of x
      xhat(n)=thetahat(n,:)*Y(:,n);
      % Update thetahat
      thetahat(n+1,:)=thetahat(n,:)+muu*Y(:,n)'*(x(n)-xhat(n));
    end
    % Shift thetahat one step so that row n corresponds to time n
    thetahat=thetahat(2:M+1,:);