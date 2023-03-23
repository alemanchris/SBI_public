%{
 The following is an example that illustrates the use of the routine:
 SBI_nrm. SBI_nrm implements the Stage Based Identification method in
 Aleman,  Busch,  Ludwig,  and  Santaeul`alia-Llopis  (2023)  
 Any comments please address them to christian.c.aleman[at]gmail.com

 EXAMPLE 2:
 Exact identification using a polynomial of degree 3 as the data generating
 process.
 This also serves as a placebo example, the policy effect is zero thus
 estimated effect should be zero.
%}
clear 
close all
clc
main_path = cd;
cd ..
sbi_path = cd;
cd(main_path);

%% Generate the Data

%-----------------------
Time = (-5:0.125:5)';
par.nt = size(Time,1);
par.tp = 33; % policy

pDC = [12,-4.5,0.3,0.2];
pDT = [7.8432,-1.2852,-0.1701,0.0729];


yC = polyval([pDC(4),pDC(3),pDC(2),pDC(1)],Time);
yT = polyval([pDT(4),pDT(3),pDT(2),pDT(1)],Time);

%-----------------------
%% Compute analitical solution

[lgth,~] = size(Time);
par.time = (1:1:lgth);
pDC_del = polyfit(par.time(1:par.tp),yC(1:par.tp),3);
pDT_del = polyfit(par.time(1:par.tp),yT(1:par.tp),3);
pDC = [pDC_del(4),pDC_del(3),pDC_del(2),pDC_del(1)];
pDT = [pDT_del(4),pDT_del(3),pDT_del(2),pDT_del(1)];
% Analitical solution 
auxh2 = (pDT(4)/pDC(4))*(pDC(2)-((1/3)*(pDC(3)^2)/pDC(4)));
auxh4 = (1/3)*(pDT(3))^2/pDT(4)-pDT(2);


psi1 = (auxh2/(-auxh4))^0.5;
psi0 = (1/3)*(pDT(3)/pDT(4)*psi1)-(1/3)*(pDC(3)/pDC(4)); 
alt_w1 = pDT(4)/(pDC(4).*psi1^3);
alt_w0 = pDT(1)-(alt_w1*(pDC(1)+pDC(2)*psi0+pDC(3)*psi0^2+pDC(4)*psi0^3));
% Get the inverse
om(1) = alt_w0;
om(2) = alt_w1;
phi(1) = (-psi0/psi1);
phi(2) = 1/psi1;

% Graph series

figure(100)
hold on
plot(Time,yT(:,1),'r-');
plot(Time,yC(:,1),'b-');
xline(Time(par.tp),'k-');
legend('RED','BLUE','Policy')

% Mapping C to T
figure(101)
hold on
plot(par.time,yT(:,1),'r-');
plot(par.time,yC(:,1),'b-');
plot(phi(1)+par.time.*phi(2),om(1)+yC(:,1).*om(2),'bx--');
xline(par.time(par.tp),'k-');
legend('RED','BLUE','Policy')

dert = -phi(1)/phi(2)+par.time./phi(2); 
% Mapping T to C
figure(102)
hold on
plot(par.time,yT(:,1),'r-');
plot(par.time,yC(:,1),'b-');
plot(-phi(1)/phi(2)+par.time./phi(2),-om(1)/om(2)+yT(:,1)./om(2),'rx--');
xline(par.time(par.tp),'k-');
legend('RED','BLUE','Policy')

%% Conduct SBI
%% Required inputs
% Asign Region Time Series
idta(:,1) = yC(:,1);
idta(:,2) = yT(:,1);
% Set Time Variable
itm = Time;
% Set policy implementation year/date
itp = Time(par.tp);
%% Optional inputs
% Choose operating system 1: Windows Default 2: Linux
iOS = 2;
% Region names
irnam = {'REST';'REG1'}; % Write REG1 not REG_1, number of characters must be equal
% Name your Outcome variable
ionam = ['My Outcome'];
% Name the units your time variable ('year','day','time')
itnam = ['Year'];
% Mapping (1:linear (Default) or 2:quadratic)
inmts = 1;
% Level Adjustment Mapping (1:Proportional or 2:(Default)Proportional + Additive)
inmlv = 2;
% Smoothing step 1: (Default)Smoothing 0: No smoothing
ismo = 0;
% Smoother
itsmo = [];
% Boostrap % 1:Boostrap(Default) 0: Skip Boostrap Step
ib = 0; 
%{
     0: Moving average, with windown = 9;
     1: Interpolation smoother CSAPS (DEFAULT) with smoothing parameter = 0.0008
     2: Chebyshev nodes with Cheby Regression
     3: B-Splines
     4: HP Filter (lambda = 50)
%}

cd(sbi_path);
SBI_nrm(idta,itm,itp,inmts,inmlv,ismo,itsmo,[],[],[],ib,[],[],irnam,itnam,ionam,[],main_path,iOS);
cd(main_path)
%% Compare mapping coeffs with analitical solution in a table

% Load results from SBI
load('output/results_table_F')
phi_SBI(1,1) = table1.estimate(7);
phi_SBI(2,1) = table1.estimate(8);
phi_SBI(3,1) = table1.estimate(5);
phi_SBI(4,1) = table1.estimate(6);
coeffs = ["Omega_0";"Omega_1";"Phi_0";"Phi_1"];
varNames = ["Coeffs","Analitical","SBI_norm"];
comp_table = table(coeffs,[om(1);om(2);phi(1);phi(2)],phi_SBI,'VariableNames',varNames);
disp(comp_table)



