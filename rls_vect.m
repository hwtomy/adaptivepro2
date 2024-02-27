function [thetahat,xhat]=rls_vect(x,Y,lambda)

    % [thetahat,xhat]=rls_vect(x,Y,lambda)
    %
    %	x			- Data sequence
    %	Y			- Data vector sequence
    %	N			- Dimension of the parameter vector
    %	lambda			- Forgetting factor
    %	thetahat		- Matrix with estimates of theta. 
    %				  Row n corresponds to time n-1
    %	xhat			- Estimate of x for n=1
    %
    %
    %
    %  rls: Recursive Least-Squares Estimation
    %
    % 	Estimator: xhat(n)=Y^{T}(:,n)thetahat(n-1)
    %
    %	thetahat is estimated using RLS. 
    %
    %	Initalization:	P(0)=10000 I, thetahat(0)=0
    %
    %     
    %     Author: 
    %
    
    % Initialize P, xhat and thetahat
    
    
    [N,M]=size(Y);
    xhat=zeros(M,1);
    thetahat=zeros(M,N);
    P=1e4*eye(N);
    
    % Loop
    
    for n=1:M,
    
      % Generate Y(n). Set elements of Y that does not exist to zero
      
    
      % Estimate of x
      xhat(n)=thetahat(n,:)*Y(:,n);
      
    
      % Update K
      K=P*Y(:,n)/(lambda+Y(:,n)'*P*Y(:,n));
      
      % Update P
      P=1/lambda*(P-K*Y(:,n)'*P);
      
      % Update the n+1 row in the matrix thetahat which in the 
      % notation in the Lecture Notes corresponds to thetahat(n)
      thetahat(n+1,:)=thetahat(n,:)+K'*(x(n)-xhat(n));
    end
    
    % Shift thetahat one step so that row n corresponds to time n
    
    thetahat=thetahat(2:M+1,:);