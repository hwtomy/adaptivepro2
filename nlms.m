function w = nlms(s,x,mu,order)
c = 1e-9;
N = length(s);
w = zeros(order, N);

for n = order:N
    s_n = s(n:-1:n-order+1);
    
    w(:,n) = w(:,n-1) + mu/(c+s_n'*s_n)*s_n*(x(n)- transpose(s_n)*w(:,n-1));   
end


