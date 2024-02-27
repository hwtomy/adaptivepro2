function w = adag(s,x,mu,order)
%%This function is the  least mean square algorithm using adagrad to update
%%learning rate. The core is introducing second moment of gradient to
%%reflect the changes.
N = length(s);
w = zeros(order, N);
r = zeros(1,N);
e = 1e-9;
G = zeros(1,order);
for n = order:N
     s_n = s(n:-1:n-order+1);
     r(n) = x(n)- transpose(s_n)*w(:,n-1);
     grad = -2*s_n'*r(n);
     G = G+grad.*grad;
     learning = mu./(sqrt(G)+e);

     w(:,n) = w(:,n-1)-(learning.*grad)';  
end

end
