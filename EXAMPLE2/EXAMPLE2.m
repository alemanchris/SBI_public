%{
 -----------
 EXAMPLE 2:
 -----------
 Standalone, just run this file.

 This example illustrates the use of the SBI_nrm.m  routine descrived in  
 Aleman,  Busch,  Ludwig,  and  Santaeul`alia-Llopis  (2023) 
 For the complete description of the SBI_nrm.m functionalities see the README file.
 

 -----------
 DESCRIPTION: 
 -----------
 Exact identification using a polynomial of degree 3 as the data generating
 process.
 This also serves as a placebo example, the policy effect is zero thus
 estimated effect should be zero.

%}
clear 
close all
clc

%% Path information:
main_path = cd;
cd ..
sbi_path = cd;
cd(main_path);

%% Generate the Data
Time = (-5:0.125:5)';
par.nt = size(Time,1);
par.tp = 33; % policy

pDC = [12,-4.5,0.3,0.2];
pDT = [7.8432,-1.2852,-0.1701,0.0729];

yC = polyval([pDC(4),pDC(3),pDC(2),pDC(1)],Time);
yT = polyval([pDT(4),pDT(3),pDT(2),pDT(1)],Time);

%% Compute analytical solution
[lgth,~] = size(Time);
par.time = (1:1:lgth);
pDC_del = polyfit(par.time(1:par.tp),yC(1:par.tp),3);
pDT_del = polyfit(par.time(1:par.tp),yT(1:par.tp),3);
pDC = [pDC_del(4),pDC_del(3),pDC_del(2),pDC_del(1)];
pDT = [pDT_del(4),pDT_del(3),pDT_del(2),pDT_del(1)];

% Analytical solution 
auxh2 = (pDT(4)/pDC(4))*(pDC(2)-((1/3)*(pDC(3)^2)/pDC(4)));
auxh4 = (1/3)*(pDT(3))^2/pDT(4)-pDT(2);
psi1 = (auxh2/(-auxh4))^0.5;
psi0 = (1/3)*(pDT(3)/pDT(4)*psi1)-(1/3)*(pDC(3)/pDC(4)); 
alt_w1 = pDT(4)/(pDC(4).*psi1^3);
alt_w0 = pDT(1)-(alt_w1*(pDC(1)+pDC(2)*psi0+pDC(3)*psi0^2+pDC(4)*psi0^3));

% Get the inverse
om(1) = alt_w0;
om(2) = alt_w1;
psi(1) = (-psi0/psi1);
psi(2) = 1/psi1;

%% Plot Analytical Solutions
plot_analytical(Time,yT,yC,par,om,psi) % Plot analytical solution
pause(10)
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
% Region names
irnam = {'REST';'REG1'}; % Write REG1 not REG_1, number of characters must be equal
% Name your Outcome variable
ionam = ['My Outcome'];
% Name the units your time variable ('year','day','time')
itnam = ['Year'];
% Mapping (1:linear (DEFAULT) or 2:quadratic)
inmts = 1;
% Level Adjustment Mapping (1:Proportional or 2:(DEFAULT)Proportional + Additive)
inmlv = 2;
% Smoothing step (0:no smoothing  or 1: (DEFAULT) apply smoother to pre_policy series)
ismo = 0;
% Smoother
itsmo = [];
% Boostrap % 1:Boostrap(DEFAULT) 0: Skip Boostrap Step
ib = 0; 
% Show graphs % 0:Dont show(DEFAULT) 1: Show benchmark graphs
ifv = 0; 
%{
% itsmo: Choose the smoother
     0: Moving average, with windown = 9;
     1: Interpolation smoother CSAPS (DEFAULT) with smoothing parameter = 0.003
     2: Polynomial Regression: Monomial or Chebyshev basis (DEFAULT Cheby basis)
     3: B-Splines
     4: HP Filter (lambda = 50)
%}

cd(sbi_path);
SBI_nrm(idta,itm,itp,inmts,inmlv,ismo,itsmo,[],[],[],ib,[],[],irnam,itnam,ionam,[],main_path,ifv);
cd(main_path)
%% Compare mapping coeffs with analytical solution in a table

disp('*******************************************')
disp('Comparison Table Analytical Solution vs SBI')
disp('*******************************************')
% Load results from SBI
load('output/results_table_F')
psi_SBI(1,1) = table1.estimate2(6);
psi_SBI(2,1) = table1.estimate2(7);
psi_SBI(3,1) = table1.estimate2(4);
psi_SBI(4,1) = table1.estimate2(5);
coeffs = ["Omega_0";"Omega_1";"psi_0";"psi_1"];
varNames = ["Coeffs","Analytical","SBI_norm"];
comp_table = table(coeffs,[om(1);om(2);psi(1);psi(2)],psi_SBI,'VariableNames',varNames);
disp(comp_table)

%% Graph Analytical Solution

% ------------------------------------------------------------------------
%
function plot_analytical(Time,yT,yC,par,om,psi)
% Graph series
figure(100)
hold on
plot(Time,yT(:,1),'r-','linewidth',1.2);
plot(Time,yC(:,1),'b-','linewidth',1.2);
xline(Time(par.tp),'k-','linewidth',1.1);
legend('RED','BLUE','Policy')

% Mapping C to T
figure(101)
hold on
plot(par.time,yT(:,1),'r-','linewidth',1.2);
plot(par.time,yC(:,1),'b-','linewidth',1.2);
plot(psi(1)+par.time.*psi(2),om(1)+yC(:,1).*om(2),'bx--');
xline(par.time(par.tp),'k-','linewidth',1.1);
title('C to T')
legend('RED','BLUE','BLUE NORM','Policy')

dert = -psi(1)/psi(2)+par.time./psi(2); 
% Mapping T to C
figure(102)
hold on
plot(par.time,yT(:,1),'r-','linewidth',1.2);
plot(par.time,yC(:,1),'b-','linewidth',1.2);
plot(-psi(1)/psi(2)+par.time./psi(2),-om(1)/om(2)+yT(:,1)./om(2),'rx--');
xline(par.time(par.tp),'k-','linewidth',1.1);
title('T to C')
legend('RED','BLUE','RED NORM','Policy')
end
%}


