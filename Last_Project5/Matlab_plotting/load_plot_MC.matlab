clear;
%------------------

%------------------
forwardEuler= 'forwardEuler.txt';
crankNicolson= 'crankNicolson.txt';
backwardEuler= 'backwardEuler.txt';
analytical='analytical.txt';
delimiterIn = ' ';
headerlinesIn = 1;
explicit = importdata(forwardEuler,delimiterIn,headerlinesIn);
implicit = importdata(backwardEuler,delimiterIn,headerlinesIn);
cn_scheme = importdata(crankNicolson,delimiterIn,headerlinesIn);
analytic = importdata(analytical,delimiterIn,headerlinesIn);
dim_explicit= size(explicit)
% u coresponding to time 0.05, 
%Note: Delta_x =0.1 yielding Delta_t=0.005
n_x1= length(explicit.data( 1,:))
n_t1= length(explicit.data( :,1))
x_exp = linspace(0,1,n_x1) ;
n_x2= length(implicit.data( 1,:))
%n_t= length(implicit.data( :,1))
x_imp = linspace(0,1,n_x2) ;
n_x3= length(cn_scheme.data( 1,:))
%n_t= length(cn_scheme.data( :,1))
x_cn = linspace(0,1,n_x3) ;
n_x4= length(analytic.data( 1,:))
%n_t= length(analytic.data( :,1))
x_ana = linspace(0,1,n_x4) ;
one=ones(1,n_x1);
prompt='Write the spacial step of the data: ';
Dx=input(prompt)
Dt=Dx*Dx


%This is for dx=0.1
nt1=int8(0.05/Dt)
nt2=int8(0.5/Dt)
%end dx=0.1

%This is for dx=0.01
%nt1=int8(0.00005/Dt)
%nt2=int8(0.9/Dt)
%end dx=0.01
%-------The prosent of deiviation from the analytical error
%This method is not the best for a visual understaning of the behavoir of
%the diffirent method. because when compared with the analytical solution
%it get affect by the change of the analytical value. so the diviation
%seems like it get bigger while in reality it get smaller and smaller.
% for n=1:n_x1
%     
%     anai=analytic.data(nt2,n);
%     ana_invs=1.0/(analytic.data(nt2,n));
%     if anai==0 | cn_scheme.data(nt2, n)==0
%         error_cn2(n)=0;
%        % prompt='Write the spacial step of the data: ';
%         %ss=input(prompt)
%     else
%         error_cn2(n)=1-(cn_scheme.data(nt2, n)*ana_invs);
%          
%     end
%     %--
%     if explicit.data(nt2, n)==0 | anai==0
%         error_exp2(n)=0;
%     else
%        % explicit_n= explicit.data(nt2, n)
%         %ana_invs
%         error_exp2(n)=1-explicit.data(nt2, n)*ana_invs;
%     end
%     %--
%     if anai==0 | implicit.data(nt2, n)==0
%         error_imp2(n)=0;
%     else
%         error_imp2(n)=1-(implicit.data(nt2, n)*ana_invs);
%     end
% end
% 
% %_------------------------
% for n=1:n_x1
%     
%     anai=analytic.data(nt1,n);
%     ana_invs=1.0/(analytic.data(nt1,n));
%     if anai==0 | cn_scheme.data(nt1, n)==0
%         error_cn1(n)=0;
%        % prompt='Write the spacial step of the data: ';
%         %ss=input(prompt)
%     else
%         error_cn1(n)=1-(cn_scheme.data(nt1, n)*ana_invs);
%     end
%     %--
%     if explicit.data(nt1, n)==0 | anai==0
%         error_exp1(n)=0;
%     else
%        % explicit_n= explicit.data(nt2, n)
%         %ana_invs
%         error_exp1(n)=1-explicit.data(nt1, n)*ana_invs;
%     end
%     %--
%     if anai==0 | implicit.data(nt1, n)==0
%         error_imp1(n)=0;
%     else
%         error_imp1(n)=1-(implicit.data(nt1, n)*ana_invs);
%     end
% end

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
%cn=cn_scheme.data(nt2, :);
%ana=analytic.data(nt2, :);
%an_1=1./analytic.data(nt2, :)
%an_pw1=analytic.data(nt2, :).^(-1)


%(explicit.data(nt2, :)./analytic.data(nt2, :))-one
%(implicit.data(nt2, :)./analytic.data(nt2, :))-one


%error_cn=(cn_scheme.data(nt2, :)./analytic.data(nt2, :))-one
% Energy=Avr.data(:, 8);
% Energy=Energy(200000:L);
% a = unique(Energy);
% Efrequency = [a,histc(Energy(:),a)]
% Efrequencyprosent= [a,histc(Energy(:),a)/length(Energy)]
% 
% 
