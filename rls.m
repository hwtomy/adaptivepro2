function w=rls(s,x,order)

N = length(s);
P = .05*eye(order);
lambda = 0.999999;
w = zeros(order, N);

for n = order:N
  Y = s(n:-1:n-order+1);

  K=P*Y/(lambda+Y'*P*Y);  
  P=1/lambda*(P-K*Y'*P);
  
  w(:,n)=w(:,n-1)+K*(x(n)-Y'*w(:,n-1));
end