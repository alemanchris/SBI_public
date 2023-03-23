%{
 EXAMPLE 1:
 Uses example data to run SBI
%}
clear 
close all
clc

main_path = cd;
cd ..
sbi_path = cd;
cd(main_path);

% Load Data
load([main_path,'/input/exampledata'])
%% Required inputs
% Asign Region Time Series
idta(:,1) = data(:,3);
idta(:,2) = data(:,1);
% Set Time Variable
itm = data_time;
% Set policy implementation year/date
itp = 1960;
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
inmlv = 1;
% Smoothing step, use default 
ismo = [];
% Smoother, use default
itsmo = [];
% Boostrap % 1:Boostrap (Default) 0: Skip Boostrap Step
ib = []; 
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


