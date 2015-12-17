clear;
%------------------

%------------------
mc_data= 'density.dat';
delimiterIn = ' ';
headerlinesIn = 1;
data = importdata(mc_data,delimiterIn,headerlinesIn);



% %This is for dx=0.1
% nt1=int8(0.05/Dt)
% nt2=int8(0.5/Dt)
% %end dx=0.1

%This is for dx=0.01
nt1=int8(0.00005/Dt)
nt2=int8(0.9/Dt)
%end dx=0.01
%-----The diviation from the analytical error
for n=1:n_x1
    anai=analytic.data(nt2,n);
    error_cn2(n)= cn_scheme.data(nt2, n)-anai;
    error_exp2(n)=explicit.data(nt2, n)-anai;
    error_imp2(n)=implicit.data(nt2, n)-anai;
end

%_------------------------
for n=1:n_x1
    anai=(analytic.data(nt1,n));
    error_cn1(n)=cn_scheme.data(nt1, n)-anai;
    error_exp1(n)=explicit.data(nt1, n)-anai;
    error_imp1(n)=(implicit.data(nt1, n)-anai);
end
%---------
figure(1)
plot(x_exp,explicit.data(nt1, :))
hold on
plot(x_imp,implicit.data(nt1, :))
hold on
plot(x_cn,cn_scheme.data(nt1, :),'--')
hold on
plot(x_ana,analytic.data(nt1, :),':')
    
title('Relative transmitters consentration over synaptic cleft')
xlabel('Synaptic cleft length ')
ylabel('Relativ consentration')
legend('Explicit Scheme','Implicit Scheme','Crank-Nicolson','Analytical Solution')
hold off
% %----------
figure(2)
plot(x_exp,explicit.data(nt2, :))
hold on
plot(x_imp,implicit.data(nt2, :))
hold on
plot(x_cn,cn_scheme.data(nt2, :),'--')
hold on
plot(x_ana,analytic.data(nt2, :),':')
    
title('Relative transmitters consentration over synaptic cleft')
xlabel('Synaptic cleft length ')
ylabel('Relativ consentration')
legend('Explicit Scheme','Implicit Scheme','Crank-Nicolson','Analytical Solution')
hold off

%Relative error:

% %---------
figure(3)
plot(x_exp,error_exp1)
hold on
plot(x_imp,error_imp1)
hold on
plot(x_cn,error_cn1,'--')
    
title('Numerical deviation of the consentration relative to the analyticiall solution')
xlabel('Synaptic cleft length ')
ylabel('Relativ error')
legend('Explicit Scheme','Implicit Scheme','Crank-Nicolson')
hold off
%----------
figure(4)
plot(x_exp,error_exp2)
hold on
plot(x_imp,error_imp2)
hold on
plot(x_cn,error_cn2,'--')

title('Numerical deviation of the consentration relative to the analyticiall solution')
xlabel('Synaptic cleft length ')
ylabel('Relativ error')
legend('Explicit Scheme','Implicit Scheme','Crank-Nicolson')
hold off   
