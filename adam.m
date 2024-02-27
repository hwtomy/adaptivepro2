function w = adam(s,x,mu,order,b1,b2)
N = length(s);
w = zeros(order, N);
r = zeros(1,N);
e = 1e-9;
g1 = zeros(1,order);
g2 = zeros(1,order);
t= 0;

for n = order:N
    s_n = s(n:-1:n-order+1);
     r(n) = x(n)- transpose(s_n)*w(:,n-1);
     grad = -2*s_n'*r(n);

     t=t+1;
     g1 = g1*b1+(1-b1)*grad;
     g2 = g2*b2+(1-b2)*grad.*grad;

     g1_h = g1/(1-b1^t);
     g2_h = g2/(1-b2^t);
     learning = mu./(sqrt(g2_h)+e);
     
     w(:,n) = w(:,n-1)-(learning.*g1_h)';  

end
end
