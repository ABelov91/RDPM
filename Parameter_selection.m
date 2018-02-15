reaction_number = 1;
[E,u,delta] = reaction_choice(reaction_number);

N = length(E); % число экспериментальных точек

Fund_period = 20; % число пар синусов-косинусов основного периода
J_max = 2*Fund_period + 4; % максимальное число аппроксимирующих слагаемых
terms_num = zeros(1,Fund_period);

E_min = min(E); E_max = max(E);
x = pi*(E - E_min)/(E_max - E_min) - pi/2;

Alpha_max = log10(0.5); Alpha_min = -10;
K = 10;
Alpha = zeros(1,K);

residual = zeros(1,(J_max-4)/2);
res  = zeros(K+1,(J_max-4)/2);
cond = zeros(K+1,(J_max-4)/2);

for k = 1:K+1
        Alpha(k) = 10^(Alpha_min + (Alpha_max - Alpha_min)*(k-1)/K);
        Beta0 = 3;
        Beta = Beta0^2*0.7;
        Reg_params = [Alpha(k), Beta, 0];  
    for J = 1:(J_max-4)/2
        
        terms_num(J) = 2*(J-1)+4;
        [fit,cond(k,J),res(k,J)] = LSQ(x,u,delta, terms_num(J), Reg_params);

    end
end

figure; hold on;
v1 = [-0.325,-0.32,-0.3]; v2 = -0.3:0.02:0;
v = cat(2,v1,v2);
[C1,h1] = contour(log10(terms_num-3),log10(Alpha),log10(res),v,'r','ShowText','on');
clabel(C1,h1,'FontSize',18,'Color','black')
v1 = 1:2:15;
[C2,h2] = contour(log10(terms_num-3),log10(Alpha),log10(cond),v1,'b','ShowText','on');
clabel(C2,h2,'FontSize',18,'Color','black')