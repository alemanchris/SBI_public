%{
 -----------
 EXAMPLE 3:
 -----------
 Standalone, just run this file.

 The following is an example of the use of the routine:
 SBI_nrm, that implements the Stage Based Identification methods in
 Aleman,  Busch,  Ludwig,  and  Santaeul`alia-Llopis  (2020)  
 For the complete description of the SBI_nrm.m functionalities see the README file.

 -----------
 DESCRIPTION: 
 -----------

 Exact identification using GLF(Generalized Logistic Func) with as the data generatig process
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

%% Generate data
% Data generating process GLF
T = 200; 
Time = (1:T)';
tp = 70;        % Placebo policy date
theta1 = [0,35,0.088,90]; 
theta2 = [0,40,0.09,77]; 
[GT,gT] = GLF(theta1,Time);
gT2 = diff(GT); % this is an approximation of the derivative.  
[GC,gC] = GLF(theta2,Time);
gC2 = diff(GC);

%% Analytical Solutions to the normalization taken from Aleman,  Busch,  Ludwig,  and  Santaeul`alia-Llopis  (2020)
an_psi = NaN(4,1);
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

an_psi(2) = (theta0T-theta1T)/(theta0C-theta1C);
an_psi(1) = theta0T-theta0C*an_psi(2);
an_psi(3) = (theta2C-theta2T)/theta3C;
an_psi(4) = theta3T/theta3C;

% Get the inverse
om(1) = an_psi(1);
om(2) = an_psi(2);
psi(1) = (-an_psi(3)/an_psi(4));
psi(2) = 1/an_psi(4);

%% Plot Analytical Solutions
plot_analytical_GLF(Time,GT,GC,tp,psi,om) % Normalization in Levels
pause(10)

%% PART 1: Run SBI for GLF
%% Required inputs
% Asign Region Time Series
idta(:,1) = GC;
idta(:,2) = GT;
% Set Time Variable
itm = Time;
% Set policy implementation year/date
itp = tp;
%% Optional inputs
% Region names
irnam = {'BLU';'RED'}; % Write REG1 not REG_1, number of characters must be equal
% Name your Outcome variable
ionam = ['y'];
% Name the units your time variable ('year','day','time')
itnam = ['time'];
% Custom Name Prefix for tables and figures
icusnam = ['GLF'];
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
SBI_nrm(idta,itm,itp,inmts,inmlv,ismo,itsmo,[],[],[],ib,[],[],irnam,itnam,ionam,icusnam,main_path,ifv);
cd(main_path);


%% Compare mapping coeffs with analytical solution in a table
disp('***************************************************')
disp('Comparison Table 1: GLF, analytical solution vs SBI')
disp('***************************************************')
% Load results from SBI
load('output/results_table_GLF')
psi_SBI(1,1) = table1.estimate2(6);
psi_SBI(2,1) = table1.estimate2(7);
psi_SBI(3,1) = table1.estimate2(4);
psi_SBI(4,1) = table1.estimate2(5);
coeffs = ["Omega_0";"Omega_1";"psi_0";"psi_1"];
varNames = ["Coeffs","Analytical","SBI_norm"];
comp_table = table(coeffs,[om(1);om(2);psi(1);psi(2)],psi_SBI,'VariableNames',varNames);
disp(comp_table)




% ------------------------------------------------------------------------
function [y,dydx] = GLF(theta,x)

y = (theta(2)-theta(1))./(1+exp(-theta(3).*( x - theta(4))))+theta(1); 
dydx = -(theta(3).*(theta(1) - theta(2)).*exp(theta(3).*(theta(4) + x)))./(exp(theta(3)*theta(4)) + exp(theta(3).*x)).^2;
end
% ------------------------------------------------------------------------
function plot_analytical_GLF(Time,GT,GC,tp,psi,om)

% Graph series
figure(100)
hold on
plot(Time,GT,'r-','linewidth',1.2);
plot(Time,GC,'b-','linewidth',1.2);
xline(tp,'k-','linewidth',1.1);
title('GLF')
legend('RED','BLUE','Policy date')

% Mapping C to T

figure(101)
hold on
plot(Time,GT(:,1),'r-','linewidth',1.2);
plot(Time,GC(:,1),'b-','linewidth',1.2);
plot(psi(1)+Time.*psi(2),om(1)+GC(:,1).*om(2),'bx--');
xline(tp,'k-','linewidth',1.1);
title('C to T')
legend('RED','BLUE','BLUE Norm','Policy')

% Mapping T to C

figure(102)
hold on
plot(Time,GT(:,1),'r-','linewidth',1.2);
plot(Time,GC(:,1),'b-','linewidth',1.2);
plot(-psi(1)/psi(2)+Time./psi(2),-om(1)/om(2)+GT(:,1)./om(2),'rx--');
xline(tp,'k-','linewidth',1.1);
title('T to C')
legend('RED','BLUE','RED Norm','Policy')

end
% ------------------------------------------------------------------------

