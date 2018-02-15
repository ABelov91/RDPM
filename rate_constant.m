function [T,K] = rate_constant(E_min0,E_max0, reaction_choice,fit)
if( reaction_choice == 1 )    
    B = 31.396;
    ma = 1.673e-24;
elseif( reaction_choice == 2 )
    B = 31.396;
    ma = 1.673e-24;
elseif( reaction_choice == 3 )
    B = 34.392;
    ma = 1.2*1.673e-24;
elseif( reaction_choice == 4 )
    B = 68.784;
    ma = 1.2*1.673e-24;
end

c = 2;
E_grid = @(x)( x/(1 - x^2)^c );
E_step = @(x)( ( 1 + x^2*(2*c-1) )/( 1 - x^2 )^(c+1) );

N_T = 100;
T = zeros(1,N_T+1);
T_min = log10(0.01);
T_max = log10(2000);
for n = 1:N_T+1
    T(n) = 10^( T_min + (T_max - T_min)*(n-1)/N_T );
end

K = zeros(1,N_T+1);

N_E = 500; % чило узлов начальной сетки

E = zeros(1,N_E);
h_E = zeros(1,N_E);
xi = zeros(1,N_E);
sigma = zeros(1,N_E);
f = zeros(1,N_E);

x = zeros(1,N_E);
E_temp = zeros(1,N_E);

for n = 1:N_T+1
    for k = 1:N_E
        xi(k) = (k-0.5)/(N_E);
        E(k) = E_grid(xi(k));
        h_E(k) = E_step(xi(k))/N_E;
        E_temp(k) = E(k);
    end

    for k = 1:N_E
        x(k) = pi*(log10(E(k)) - E_min0)/(E_max0 - E_min0) - pi/2;
    end

    J = length(fit);

    phi0 = zeros(J,N_E);
    phi0(1,:) = sin(x);
    phi0(2,:) = cos(x);
    phi0(3,:) = sin(3*x);
    phi0(4,:) = ones(1,N_E);
    for j = 3:J/2
        phi0(2*j-1,:) = sin( (2*j-4)*x );
        phi0( 2*j ,:) = cos( (2*j-4)*x );
    end

    for k = 1:N_E-1
        E01 = log10(E(k));
        E02 = log10(E(k+1));
        if( ( E01 - E_min0 )*( E02 - E_min0 ) < 0 )
            stitch_num1 = k;
        end
        if( ( E_max0 - E01 )*( E_max0 - E02 ) < 0 )
            stitch_num2 = k;
        end
    end
    x_temp1 = x(stitch_num1);
    x_temp2 = x(stitch_num2);
    E_temp1 = log10( E(stitch_num1) );
    E_temp2 = log10( E(stitch_num2) );

    value1 = 0;
    for j = 1:J
        value1 = value1 + fit(j)*phi0(j,stitch_num1);
    end
    slope1 = fit(1)*cos(x_temp1) - fit(2)*sin(x_temp1) + 3*fit(3)*cos(3*x_temp1);
    for j = 3:J/2
        slope1 = slope1 + (2*j-4)*fit(2*j-1)*cos( (2*j-4)*x_temp1 ) - (2*j-4)*fit(2*j)*sin( (2*j-4)*x_temp1 );
    end
    value2 = 0;
    for j = 1:J
        value2 = value2 + fit(j)*phi0(j,stitch_num2);
    end
    slope2 = fit(1)*cos(x_temp2) - fit(2)*sin(x_temp2) + 3*fit(3)*cos(3*x_temp2);
    for j = 3:J/2
        slope2 = slope2 + (2*j-4)*fit(2*j-1)*cos( (2*j-4)*x_temp2 ) - (2*j-4)*fit(2*j)*sin( (2*j-4)*x_temp2 );
    end

    approximation0 = zeros(1,N_E);
    for k = 1:N_E
        E0 = log10(E(k));
        if( E0 < E_temp1 )
            approximation0(k) = value1 + slope1*(E0 - E_temp1);
        end
        if( E0 >= E_temp1 && E0 <= E_temp2 )
            for j = 1:J
                approximation0(k) = approximation0(k) + fit(j)*phi0(j,k);
            end
        end
        if( E0 > E_temp2 )
            approximation0(k) = value2 + slope2*(E0 - E_temp2);
        end
    end
    approximation0 = 10.^approximation0/10^3;

    for k = 1:N_E
        sigma(k) = approximation0(k)/E(k)*exp( -B/sqrt( E(k) ) );            
        f(k) = 1/(T(n))^1.5*exp( -E(k)/T(n) )*E(k)*sigma(k);
    end        

    for k = 1:N_E
        K(n) = K(n) + f(k)*h_E(k);
    end
end

K = K*2e-24*sqrt( 2*1.602e-9/pi/ma );
end