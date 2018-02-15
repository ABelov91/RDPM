function [fit,conditioning,residual] = LSQ(x,u,delta, J,Alpha)
N = length(x);
phi = zeros(J,N);
phi(1,:) = sin(x);
phi(2,:) = cos(x);
phi(3,:) = sin(3*x);
phi(4,:) = ones(1,N);
for j = 3:J/2
    phi(2*j-1,:) = sin( (2*j-4)*x );
    phi( 2*j ,:) = cos( (2*j-4)*x );
end

[Gram_matrix,right_hand] = Least_squares_equations(J, N, u, phi, delta, Alpha);
conditioning = cond(Gram_matrix);
fit = Gram_matrix\right_hand';

approximation = zeros(1,N);
for j = 1:J
    approximation = approximation + fit(j)*phi(j,:);
end
errors = u - approximation;
residual = 0;
for j = 1:N
    residual = residual + (errors(j)./delta(j)).^2;
end
residual = sqrt( residual/(N-J) );
end