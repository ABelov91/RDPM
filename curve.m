function A = curve(x,fit)
N = length(x);
J = length(fit);

phi = zeros(J,N);
phi(1,:) = sin(x);
phi(2,:) = cos(x);
phi(3,:) = sin(3*x);
phi(4,:) = ones(1,N);
for j = 3:J/2
    phi(2*j-1,:) = sin( (2*j-4)*x );
    phi( 2*j ,:) = cos( (2*j-4)*x );
end

A = zeros(1,N);
for j = 1:J
    A = A + fit(j)*phi(j,:);
end
end