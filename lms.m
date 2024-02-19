function w = lms(s,x,mu,order)
N = length(s);
w = zeros(order, N);

for n = order:N
    s_n = s(n:-1:n-order+1);
    
    w(:,n) = w(:,n-1) + mu*s_n*(x(n)- transpose(s_n)*w(:,n-1));   
end


