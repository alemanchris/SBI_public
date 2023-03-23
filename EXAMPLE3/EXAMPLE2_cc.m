%{
 The following is an example of the use of the routine:
 SBI_nrm, that implements the Stage Based Identification methods in
 Aleman,  Busch,  Ludwig,  and  Santaeul`alia-Llopis  (2020)  
 Any comments please address them to christian.c.aleman[at]gmail.com

 EXAMPLE 2:
 Exact identification using GLF(Generalized Logistic Func) with as the data generatig process
 This also serves as a placebo example, the policy effect is zero thus
 estimated effect should be zero.
 - This example shows the perfomance of SBI on the Generaliwed logictic
 function(GLF) and its derivative(dydxGFL)
%}
clear 
close all
clc

%% path information:
main_path = cd;
cd ..
sbi_path = cd;
cd(main_path);

%% Generate data
% Data generating process GLF
T = 200; 
Time = (1:T)';
tp = 70;        % Placebo policy date
theta1 = [0,35,0.088,90]; 
theta2 = [0,40,0.09,77]; 
GT = GLF(theta1,Time);
gT = diff(GT);
GC = GLF(theta2,Time);
gC = diff(GC);

%% Analitical Solutions to the normalization taken from Aleman,  Busch,  Ludwig,  and  Santaeul`alia-Llopis  (2020)
an_phi = NaN(4,1);
theta1(4) = theta1(3)*theta1(4);
theta2(4) = theta2(3)*theta2(4);

theta0C = theta2(1);
theta1C = theta2(2);
theta2C = theta2(4);
theta3C = theta2(3);

theta0T = theta1(1);
theta1T = theta1(2);
theta2T = theta1(4);
theta3T = theta1(3);

an_phi(2) = (theta0T-theta1T)/(theta0C-theta1C);
an_phi(1) = theta0T-theta0C*an_phi(2);
an_phi(3) = (theta2C-theta2T)/theta3C;
an_phi(4) = theta3T/theta3C;

% Get the inverse
om(1) = an_phi(1);
om(2) = an_phi(2);
phi(1) = (-an_phi(3)/an_phi(4));
phi(2) = 1/an_phi(4);

%% Normalization in Levels
% Graph series
figure(100)
hold on
plot(Time,GT,'r-');
plot(Time,GC,'b-');
xline(tp,'k-');
title('GLF')
legend('RED','BLUE','Policy date')

% Mapping C to T

figure(101)
hold on
plot(Time,GT(:,1),'r-');
plot(Time,GC(:,1),'b-');
plot(phi(1)+Time.*phi(2),om(1)+GC(:,1).*om(2),'bx--');
xline(tp,'k-');
legend('RED','BLUE','BLUE Norm','Policy')

% Mapping T to C

figure(102)
hold on
plot(Time,GT(:,1),'r-');
plot(Time,GC(:,1),'b-');
plot(-phi(1)/phi(2)+Time./phi(2),-om(1)/om(2)+GT(:,1)./om(2),'rx--');
xline(tp,'k-');
legend('RED','BLUE','RED Norm','Policy')

%% Run SBI for GLF
% Required inputs
% Asign Region Time Series
idta(:,1) = GC;
idta(:,2) = GT;
% Set Time Variable
itm = Time;
% Set policy implementation year/date
itp = tp;
%% Optional inputs
% Region names
iOS = 2; % 1: Windows 2: Linux
% Region names
irnam = {'BLU';'RED'}; % Write REG1 not REG_1, number of characters must be equal
% Name your Outcome variable
ionam = ['y'];
% Name the units your time variable ('year','day','time')
itnam = ['time'];
% Custom Name Prefix for tables and figures
icusnam = ['GLF'];
% Mapping (1:linear (Default) or 2:quadratic)
inmts = 1;
% Level Adjustment Mapping (1:Proportional or 2:(Default)Proportional + Additive)
inmlv = 2;
% Smoothing step (0:no smoothing or 1:apply smoother to series)
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
SBI_nrm(idta,itm,itp,inmts,inmlv,ismo,itsmo,[],[],[],ib,[],[],irnam,itnam,ionam,icusnam,main_path,iOS);
cd(main_path);


%% Compare mapping coeffs with analitical solution in a table

% Load results from SBI
load('output/results_table_GLF')
phi_SBI(1,1) = table1.estimate(7);
phi_SBI(2,1) = table1.estimate(8);
phi_SBI(3,1) = table1.estimate(5);
phi_SBI(4,1) = table1.estimate(6);
coeffs = ["Omega_0";"Omega_1";"Phi_0";"Phi_1"];
varNames = ["Coeffs","Analitical","SBI_norm"];
comp_table = table(coeffs,[om(1);om(2);phi(1);phi(2)],phi_SBI,'VariableNames',varNames);
disp(comp_table)

%% Normalization of the Derivative dy/dxGLF

% Plot the Series 
figure(200)
hold on
plot(Time(1:end-1),gT,'r-');
plot(Time(1:end-1),gC,'b-');
xline(tp,'k-');
title('dydx GLF')
legend('RED','BLUE','Policy date')

% Plot Analytical Solution

% Mapping C to T
figure(201)
hold on
plot(Time(1:end-1),gT(:,1),'r-');
plot(Time(1:end-1),gC(:,1),'b-');
plot(phi(1)+Time(1:end-1).*phi(2),om(1)+gC(:,1).*om(2),'bx--');
xline(tp,'k-');
legend('RED','BLUE','BLUE Norm','Policy')

% Mapping C to T
figure(202)
hold on
plot(Time(1:end-1),gT(:,1),'r-');
plot(Time(1:end-1),gC(:,1),'b-');
plot(-phi(1)/phi(2)+Time(1:end-1)./phi(2),-om(1)/om(2)+gT(:,1)./om(2),'rx--');
xline(tp,'k-');
legend('RED','BLUE','RED Norm','Policy')

%% Run SBI for dy/dxGLF
% Required inputs
% Asign Region Time Series
clear idta 
idta(:,1) = gC;
idta(:,2) = gT;
% Set Time Variable
itm = Time(1:end-1);
% Set policy implementation year/date
itp = tp;
%% Optional inputs
% Custom Name Prefix for tables and figures
icusnam = ['dydxGLF'];
%{
     0: Moving average, with windown = 9;
     1: Interpolation smoother CSAPS (DEFAULT) with smoothing parameter = 0.0008
     2: Chebyshev nodes with Cheby Regression
     3: B-Splines
     4: HP Filter (lambda = 50)
%}
%SBI_nrm(idta,itm,itp,inmts,inmlv,ismo,itsmo,[],ib,irnam,itnam,ionam);

cd(sbi_path);
SBI_nrm(idta,itm,itp,inmts,inmlv,ismo,itsmo,[],[],[],ib,[],[],irnam,itnam,ionam,icusnam,main_path,iOS);
cd(main_path);


%% Compare mapping coeffs with analitical solution in a table

% Load results from SBI
load('output/results_table_dydxGLF')
phi_SBI(1,1) = table1.estimate(7);
phi_SBI(2,1) = table1.estimate(8);
phi_SBI(3,1) = table1.estimate(5);
phi_SBI(4,1) = table1.estimate(6);
coeffs = ["Omega_0";"Omega_1";"Phi_0";"Phi_1"];
varNames = ["Coeffs","Analitical","SBI_norm"];
comp_table = table(coeffs,[om(1);om(2);phi(1);phi(2)],phi_SBI,'VariableNames',varNames);
disp(comp_table)

% ------------------------------------------------------------------------
function [y] = GLF(theta,x)

y = (theta(2)-theta(1))./(1+exp(-theta(3).*( x - theta(4))))+theta(1); 

end






