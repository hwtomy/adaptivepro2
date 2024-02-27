function w=rls(s,x,order, lambda)
%%This function is the standard recursive least mean square algorithm
N = length(s);
P = .05*eye(order);
w = zeros(order, N);

for n = order:N
  Y = s(n:-1:n-order+1);

  K=P*Y/(lambda+Y'*P*Y);  
  P=1/lambda*(P-K*Y'*P);

  w(:,n)=w(:,n-1)+K*(x(n)-Y'*w(:,n-1));
end

end