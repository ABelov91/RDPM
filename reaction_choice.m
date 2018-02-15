function [E,u,delta, J, Alpha] = reaction_choice(reaction_number)
if( reaction_number == 1 )
    file = cellstr( char( 'SF_D+D=p+T_Arnold.txt',...
                      'SF_D+D=p+T_Booth.txt',...
                      'SF_D+D=p+T_Bretscher.txt',...
                      'SF_D+D=p+T_Brolley.txt',...
                      'SF_D+D=p+T_Brown.txt',...
                      'SF_D+D=p+T_Chinese.txt',...
                      'SF_D+D=p+T_Cook.txt',...
                      'SF_D+D=p+T_Davenport.txt',...
                      'SF_D+D=p+T_Ganeev.txt',...
                      'SF_D+D=p+T_Graves.txt',...
                      'SF_D+D=p+T_Greife_1.txt',...
                      'SF_D+D=p+T_Greife_2.txt',...
                      'SF_D+D=p+T_Gruebler_1972.txt',...
                      'SF_D+D=p+T_Gruebler_1981.txt',...
                      'SF_D+D=p+T_Krauss.txt',...
                      'SF_D+D=p+T_Krauss_large_energies.txt',...
                      'SF_D+D=p+T_Leonard.txt',...
                      'SF_D+D=p+T_Mcneill.txt',...
                      'SF_D+D=p+T_Moffatt.txt',...
                      'SF_D+D=p+T_Pizzone.txt',...
                      'SF_D+D=p+T_Preston.txt',...
                      'SF_D+D=p+T_Raiola_2002.txt',...
                      'SF_D+D=p+T_Raiola_2002_1.txt',...
                      'SF_D+D=p+T_Raiola_2002_2.txt',...
                      'SF_D+D=p+T_Raiola_2004.txt',...
                      'SF_D+D=p+T_Raiola_2004_1.txt',...
                      'SF_D+D=p+T_Raiola_2005.txt',...
                      'SF_D+D=p+T_Raiola_2005_1.txt',...
                      'SF_D+D=p+T_Sanders.txt',...
                      'SF_D+D=p+T_Tumino_1.txt',...
                      'SF_D+D=p+T_Tumino_2.txt',...
                      'SF_D+D=p+T_Volkov.txt',...
                      'SF_D+D=p+T_Von-Engel.txt',...
                      'SF_D+D=p+T_Wang_Tie-Shan.txt',...
                      'SF_D+D=p+T_Wenzel.txt',...
                      'SF_D+D=p+T_Ying.txt') );
    
    J = 104;
    Alpha = [0.3,1000,0];
end
if( reaction_number == 2 )
    file = cellstr( char( 'SF_D+D=n+He3_Arnold.txt',...
                      'SF_D+D=n+He3_Belov.txt',...
                      'SF_D+D=n+He3_Booth.txt',...
                      'SF_D+D=n+He3_Brown.txt', ...
                      'SF_D+D=n+He3_Bystritsky_2008.txt',...
                      'SF_D+D=n+He3_Bystritsky_2010.txt',...
                      'SF_D+D=n+He3_Bystritsky_2012_1.txt',...
                      'SF_D+D=n+He3_Bystritsky_2012_2.txt',...
                      'SF_D+D=n+He3_Bystritsky_2012_3.txt',...
                      'SF_D+D=n+He3_Chagnon.txt',...
                      'SF_D+D=n+He3_Chinese.txt',...
                      'SF_D+D=n+He3_Davidenko.txt',...
                      'SF_D+D=n+He3_Ganeev_1.txt',...
                      'SF_D+D=n+He3_Ganeev_2.txt',...
                      'SF_D+D=n+He3_Goldberg.txt',...
                      'SF_D+D=n+He3_Greife.txt',...
                      'SF_D+D=n+He3_Hofstee.txt',...
                      'SF_D+D=n+He3_Hunter.txt',...
                      'SF_D+D=n+He3_Krauss_large_energies.txt',...
                      'SF_D+D=n+He3_Krauss.txt',...
                      'SF_D+D=n+He3_Leonard.txt',...
                      'SF_D+D=n+He3_Manley.txt',...
                      'SF_D+D=n+He3_Mcneill.txt',...
                      'SF_D+D=n+He3_Preston.txt',...
                      'SF_D+D=n+He3_Schulte.txt',...
                      'SF_D+D=n+He3_Thornton.txt',...
                      'SF_D+D=n+He3_Tumino.txt',...
                      'SF_D+D=n+He3_Ying.txt') );
                  
	J = 84;
    Alpha = [0.3,1000,0];
end
if( reaction_number == 3 )
    file = cellstr( char( 'SF_D+T=n+He4_Allan.txt',...
                          'SF_D+T=n+He4_Argo_T_beam.txt',...
                          'SF_D+T=n+He4_Arnold_full.txt',...
                          'SF_D+T=n+He4_Arnold.txt',...
                          'SF_D+T=n+He4_Balabanov_D_beam.txt',...
                          'SF_D+T=n+He4_Balabanov_T_beam.txt',...
                          'SF_D+T=n+He4_Bame.txt', ...
                          'SF_D+T=n+He4_Bretscher_T_Beam.txt',...
                          'SF_D+T=n+He4_Brown.txt', ...
                          'SF_D+T=n+He4_Chinese.txt',...
                          'SF_D+T=n+He4_Davidenko_1.txt',...
                          'SF_D+T=n+He4_Davidenko_2.txt',...
                          'SF_D+T=n+He4_Drosg_new_points.txt',...
                          'SF_D+T=n+He4_Galonsky.txt',...
                          'SF_D+T=n+He4_Goldberg.txt',...
                          'SF_D+T=n+He4_Jarmie_T_beam.txt',...
                          'SF_D+T=n+He4_Jarvis.txt',...
                          'SF_D+T=n+He4_Magiera.txt',...
                          'SF_D+T=n+He4_Stewart.txt' ) );
    
    J = 34;
    Alpha = [0.005,100,0];
end

if( reaction_number == 4 )
    file = cellstr( char( 'SF_D+He3=p+He4_Aliotta.txt',...
                      'SF_D+He3=p+He4_Arnold.txt',...
                      'SF_D+He3=p+He4_Bonner.txt', ...
                      'SF_D+He3=p+He4_Carlton_He3_beam.txt',...
                      'SF_D+He3=p+He4_Cognata.txt',...
                      'SF_D+He3=p+He4_Constantini.txt',...
                      'SF_D+He3=p+He4_Engstler_D1+ions.txt',...
                      'SF_D+He3=p+He4_Engstler_D2+ions.txt',...
                      'SF_D+He3=p+He4_Engstler_D3+ions.txt',...
                      'SF_D+He3=p+He4_Engstler_He3+ions.txt',...
                      'SF_D+He3=p+He4_Erramli.txt',...
                      'SF_D+He3=p+He4_Freimer.txt',...
                      'SF_D+He3=p+He4_Geist_1.txt',...
                      'SF_D+He3=p+He4_Geist_2.txt',...
                      'SF_D+He3=p+He4_Geist_3.txt',...
                      'SF_D+He3=p+He4_Gruebler.txt',...
                      'SF_D+He3=p+He4_Jarvis.txt',...
                      'SF_D+He3=p+He4_Krauss.txt',...
                      'SF_D+He3=p+He4_Krauss_large_energies.txt',...
                      'SF_D+He3=p+He4_Kunz_1_He3_beam.txt',...
                      'SF_D+He3=p+He4_Kunz_2_He3_beam.txt',...
                      'SF_D+He3=p+He4_Moller_He3_beam.txt',...
                      'SF_D+He3=p+He4_Schroeder.txt',...
                      'SF_D+He3=p+He4_Schroeder_1989.txt',...
                      'SF_D+He3=p+He4_Stewart.txt',...
                      'SF_D+He3=p+He4_Zhichang.txt' ) );
    
	J = 54;
    Alpha = [0.05,100,0];
end

[E,u,delta] = file_read(file);
end