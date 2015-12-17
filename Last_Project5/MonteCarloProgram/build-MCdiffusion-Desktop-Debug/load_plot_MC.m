clear;
%------------------
%------------------

%fig=openfig('matlab_plot/density_time_steps5.fig','visible')
mc_density= 'density001_timesteps5000.dat';
mc_densityG= 'gausdensity001_timesteps5000.dat';
delimiterIn = ' ';
headerlinesIn = 0;
density = importdata(mc_density,delimiterIn,headerlinesIn);
densityG = importdata(mc_densityG,delimiterIn,headerlinesIn);
density2=density(:,2)/1000
density1=density(:,1)/10
density2G=densityG(:,2)/1000
density1G=densityG(:,1)/10
figure(1)
plot(density1,density2,'--')


hold on
%figure(2)
plot(density1G,density2G)
title('Monte Carlo simulation of the diffusion')
xlabel('Synaptic cleft length ')
ylabel('Relativ density')
legend('Metropolis','Gaussian random walk')
hold off
