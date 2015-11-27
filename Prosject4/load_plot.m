filename = 'outputfilebfunctionof_mcs_1011';%the file should contain the following data
%mcs,T, <E>, C_V,<M>,chi, <|M|>, flipp

%T temperature
%E energy of the state
%M magnetizisation
%This is how it's calculated in the C++ code
 %average[0] += <E>;    average[1] += E*E;
 %average[2] += M;    average[3] += M*M; average[4] += fabs(M);


%This is the example of how do file start:
%0.0000000  -2.0000000       -NAN      1.0000000     -NAN      1.0000000
%1.0000000  -1.9956000    0.027456000 0.99890000  0.001716000  0.99890000
%2.0000000  -1.7830000    0.4393000  0.92470000   0.23348200   0.92470000
%3.0000000   -0.93000000   1.1927111 0.26560000    6.4894187    0.28440000
delimiterIn = ' ';
headerlinesIn = 1;
Avr = importdata(filename,delimiterIn,headerlinesIn);

% prompt = 'What is start temperature? ';
% T0 = input(prompt)
% prompt = 'What is end temperature? ';
% Te = input(prompt)
% prompt = 'What is the increase size(steps)? ';
% dT = input(prompt)
% T=T0:dT:Te;
%seeing what we got, the following code is not needed
%test start

L = length(Avr.data(:, 1))
%----------------------------
figure(1)
plot(Avr.data(:, 1),Avr.data(:, 3),Avr.data(:, 1),1.995982)
title('Average energy <E> as function of mcs from numerical analytic calc ')
xlabel('Antall Monte Carlo circler(mcs) ')
%K boltzman constant, T temp in K, J coupling strength')
ylabel(' <E> (J:couplings strength)')
%--------------------------
% %----------------------------
% figure(2)
% plot(Avr.data(:, 1),Avr.data(:, 4),Avr.data(:, 1),0.03208)
% title('Spesific heat capasity Cv as function of mcs, numerical analytic')
% xlabel('Antall Monte Carlo circler(mcs) ')
% %K boltzman constant, T temp in K, J coupling strength')
% ylabel('Cv (1/kT^2)')
% %---------------------------
% 
% %----------------------------
% figure(3)
% plot(Avr.data(:, 1),Avr.data(:, 5))
% title('Average magnetisation <M> as function of mcs, numerically calculated')
% xlabel('Antall Monte Carlo circler(mcs)')
% %K boltzman constant, T temp in K, J coupling strength')
% ylabel('<M> ')
% %--------------------------
% figure(4)
% plot(Avr.data(:, 1),Avr.data(:, 7),0.99866)
% title('Average abs magnetisation <|M1> as function of mcs, numerical analytic.')
% xlabel('Antall Monte Carlo circler(mcs)')
% %K boltzman constant, T temp in K, J coupling strength')
% ylabel('<|M|>')
% %----------------------------
% figure(5)
% plot(Avr.data(:, 1),Avr.data(:, 6),0.0040107 )
% title('Susceptibility {\chi} as function of mcs, numerical analytic.')
% xlabel(' Antall Monte Carlo circler(mcs)')
% %K boltzman constant, T temp in K, J coupling strength')
% ylabel('{\chi} ({\beta})')
%--------------------------
% %Function of Temperature:
% %----------------------------
% figure(1)
% plot(Avr.data(:, 2),Avr.data(:, 3))
% title('<E>')
% xlabel('Temperature in (kT/J)')
% %K boltzman constant, T temp in K, J coupling strength')
% ylabel('x')
% %--------------------------
% %----------------------------
% figure(2)
% plot(Avr.data(:, 2),Avr.data(:, 4))
% title('Cv')
% xlabel('Temperature in (kT/J)')
% %K boltzman constant, T temp in K, J coupling strength')
% ylabel('x')
% %---------------------------
% 
% %T, <E>, <E^2>,,<M>, <M^2>, <|M|>
% %----------------------------
% figure(3)
% plot(Avr.data(:, 2),Avr.data(:, 5))
% title('<M>')
% xlabel('Temp')
% %K boltzman constant, T temp in K, J coupling strength')
% ylabel('<M>')
% %--------------------------
% figure(4)
% plot(Avr.data(:, 2),Avr.data(:, 5))
% title('<|M|>')
% xlabel('Temp')
% %K boltzman constant, T temp in K, J coupling strength')
% ylabel('<|M|>')
% %----------------------------
% figure(5)
% plot(Avr.data(:, 2),Avr.data(:, 6) )
% title('Chi')
% xlabel('Temperature in (kT/J)')
% %K boltzman constant, T temp in K, J coupling strength')
% ylabel('Chi')
% %--------------------------


%End test