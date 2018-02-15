reaction_number = 1;
if( reaction_number == 1 )
    reaction_name = 'D+D=p+T';
end
if( reaction_number == 2 )
    reaction_name = 'D+D=n+He3';
end
if( reaction_number == 3 )
    reaction_name = 'D+T=n+He4';
end
if( reaction_number == 4 )
    reaction_name = 'D+He3=p+He4';
end

[E,u,delta, J, Alpha] = reaction_choice(reaction_number);

E_min = min(E); E_max = max(E);
x = pi*(E - E_min)/(E_max - E_min) - pi/2;

N = length(x);
[fit,conditioning] = LSQ(x,u,delta, J,Alpha);

x_sort = sort(x); E_sort = sort(E);

S = curve(x,fit);
S_sort = curve(x_sort,fit);

[T,K] = rate_constant(E_min,E_max, reaction_number, fit);

P = 30;

M = length(T);
A1 = zeros(P,N);
K1 = zeros(P,M);
fit1 = zeros(P,J);
for p = 1:P
    u1 = zeros(1,N);
    for n = 1:N
        u1(n) = delta(n)*randn(1) + S(n);
    end
    [fit1(p,:),conditioning1] = LSQ(x,u1,delta, J,Alpha);
    A1(p,:) = curve(x_sort,fit1(p,:));
    [~,K1(p,:)] = rate_constant(E_min,E_max, reaction_number, fit1(p,:));
end

S_mn = sum(A1,1)/P;
S_sq = sum(A1.^2,1)/P;
disp_S = S_sq - S_mn.^2;
disp_S_percent = sqrt(disp_S)./S_sort*100;

K_mn = sum(K1,1)/P;
K_sq = sum(K1.^2,1)/P;
disp_K = K_sq - K_mn.^2;
disp_K_percent = sqrt(disp_K)./K*100;

figure; hold on;
xlabel('lg E, keV'); ylabel('lg S, keV mbn')
plot(E,u,'ok','MarkerFaceColor','k','MarkerSize',3)
plot(E_sort,S_sort,'-g','LineWidth',2)
axis([0.95*min(E) 1.05*max(E) -Inf Inf])
output_file = strcat( 'SFactor_',reaction_name );
print(output_file,'-dpng','-r150')

figure;
plot(E_sort,disp_S_percent,'-k','LineWidth',1.5)
xlabel('lg E, keV')
ylabel('\delta S, %')
output_file = strcat( 'SFactor_error_',reaction_name );
print(output_file,'-dpng','-r150')

figure; hold on;
plot(log10(T),log10(K),'-k','LineWidth',1.5)
xlabel('lg T, keV')
ylabel('lg K, cm^3 s^{-1}')
output_file = strcat( 'Reaction_rate_',reaction_name );
print(output_file,'-dpng','-r150')

figure;
plot(log10(T),disp_K_percent,'-k','LineWidth',1.5)
xlabel('lg T, keV')
ylabel('\delta K, %')
output_file = strcat( 'Reaction_rate_error_',reaction_name );
print(output_file,'-dpng','-r150')