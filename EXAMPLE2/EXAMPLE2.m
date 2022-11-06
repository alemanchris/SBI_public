%{
 EXAMPLE 2:
 Exact identification using GLF as the data generatig process
 This also serves as a placebo example, the policy effect is zero thus
 estimated effect should be zero.
%}
clear 
close all
clc

%% Generate data
% Data generating process GLF
T = 200; 
Time = (1:T)';
tp = 70;        % Placebo policy date
theta1 = [35,0.088,90]; 
theta2 = [40,0.09,77]; 
GT = GLF(theta1,Time);
gT = diff(GT);
GC = GLF(theta2,Time);
gC = diff(GC);
% Graph series
figure(100)
hold on
plot(Time(1:end-1),gT,'r-');
plot(Time(1:end-1),gC,'b-');
xline(tp,'k-');
legend('RED','BLUE','Policy')
%% Analitical Solutions to the normalization taken from Aleman,  Busch,  Ludwig,  and  Santaeul`alia-Llopis  (2020)

par.beta0(1) = theta2(1);
par.beta0(2) = theta1(1);
par.beta1(1) = theta2(2);
par.beta1(2) = theta1(2);
par.beta2(1) = theta2(3);
par.beta2(2) = theta1(3);

an_phi(1,1) = (par.beta0(2)*par.beta1(2))/(par.beta0(1)*par.beta1(1));
an_phi(3,1) = par.beta1(1)/par.beta1(2);
an_phi(2,1) = par.beta2(2)-an_phi(3,1)*par.beta2(1);

%% Run SBI
% Required inputs
% Asign Region Time Series
idta(:,1) = gT;
idta(:,2) = gC;
% Set Time Variable
itm = Time(1:end-1);
% Set policy implementation year/date
itp = tp;
%% Optional inputs
% Region names
irnam = ['RED';'BLU']; % Write REG1 not REG_1, number of characters must be equal
% Name your Outcome variable
ionam = ['y'];
% Name the units your time variable ('year','day','time')
itnam = ['time'];
% Mapping (2:linear or 3:quadratic)
inm = 2;
% Smoothing step (0:no smoothing or 1:apply smoother to series)
ismo = 0; 
% Smoother
itsmo = [];
%{
     0: Moving average, with windown = 9;
     1: Interpolation smoother CSAPS (DEFAULT) with smoothing parameter = 0.0008
     2: Chebyshev nodes with Cheby Regression
     3: B-Splines
     4: HP Filter (lambda = 50)
%}
SBI_nrm(idta,itm,itp,inm,ismo,itsmo,[],[],irnam,itnam,ionam);

%% Compare mapping coeffs with analitical solution in a table

% Load results from SBI
load('output/results_table')
phi_SBI(1,1) = table1.estimate(5);
phi_SBI(2,1) = table1.estimate(6);
phi_SBI(3,1) = table1.estimate(7);
coeffs = ["Phi0";"Phi1";"Phi2"];
varNames = ["Coeffs","Analitical","SBI_norm"];
comp_table = table(coeffs,an_phi,phi_SBI,'VariableNames',varNames);
disp(comp_table)


% ------------------------------------------------------------------------
function [y] = GLF(theta,x)

y = theta(1)./(1+exp(-theta(2).*( x - theta(3)))); 

end






