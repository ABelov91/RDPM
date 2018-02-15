function sum = scalar_multiplication( u, v, delta )
sum = 0;
N = length(u);
for n = 1:N-1
    sum = sum + u(n)*v(n)/delta(n)^2;
end
end