function [Gram_matrix,right_hand] = Least_squares_equations(J, N, u, phi, delta, Alpha)
Gram_matrix = zeros(J);
right_hand = zeros(1,J);
for i = 1:J
    for j = 1:J
        Gram_matrix(i,j) = scalar_multiplication( phi(i,:), phi(j,:), delta );
    end
    right_hand(i) = scalar_multiplication( u, phi(i,:), delta );
end
Gram_matrix = Gram_matrix/(N-J-3);
right_hand = right_hand/(N-J-3);

I1 = zeros(J);
I1(1,1) = pi/2 * 1^2 * 1^2 * 2;
I1(2,2) = pi/2 * 1^2 * 1^2 * 2;
I1(3,3) = pi/2 * 3^2 * 3^2 * 2;
for k = 3:J/2
    I1(2*k-1,1) = (-1)^(k+1)*(2*k-8)/(2*k-3)/(2*k-5) * 1^2 * (2*k-4)^2 * 2;
    I1(2*k,2) = (-1)^(k+1)*2/(2*k-3)/(2*k-5) * 1^2 * (2*k-4)^2 * 2;
    I1(2*k-1,3) = (-1)^(k)*(4*k-8)/(2*k-1)/(2*k-7) * 3^2 * (2*k-4)^2 * 2;
    
    I1(1,2*k-1) = (-1)^(k+1)*(2*k-8)/(2*k-3)/(2*k-5) * (2*k-4)^2 * 1^2 * 2;
    I1(3,2*k-1) = (-1)^(k)*(4*k-8)/(2*k-1)/(2*k-7) * (2*k-4)^2 * 3^2 * 2;
    I1(2*k-1,2*k-1) = pi/2 * (2*k-4)^2 * (2*k-4)^2 * 2;
    
    I1(2,2*k) = (-1)^(k+1)*2/(2*k-3)/(2*k-5) * (2*k-4)^2 * 1^2 * 2;
    I1(2*k,2*k) = pi/2 * (2*k-4)^2 * (2*k-4)^2 * 2;
end
I1 = I1/2;

I2 = zeros(J);
I2(2,2) = 2;
for j = 3:J/2
    I2(2,2*j-1) = 2*(2*j-4)*(-1)^j;
    I2(2*j-1,2) = 2*(2*j-4)*(-1)^j;
    
    for i = 3:J/2
        I2(2*i-1,2*j-1) = 2*(2*j-4)*(2*i-4)*(-1)^(i+j);
    end
end
I2 = I2/2;

I3 = zeros(J);
I3(1,1) = 2;
I3(3,3) = 162;
I3(1,3) = -18;
I3(3,1) = -18;
for j = 3:J/2
    I3(1,2*j) = -2*(2*j-4)^2*(-1)^j;
    I3(3,2*j) = 18*(2*j-4)^2*(-1)^j;
    
    I3(2*j,1) = -2*(2*j-4)^2*(-1)^j;
    I3(2*j,3) = 18*(2*j-4)^2*(-1)^j;
    for i = 3:J/2
        I3(2*j,2*i) = 2*(2*j-4)^2*(2*i-4)^2*(-1)^(i+j);
    end
end
I3 = I3/2;

Gram_matrix = Gram_matrix + Alpha(1)*I1 + Alpha(2)*I2 + Alpha(3)*I3;
end