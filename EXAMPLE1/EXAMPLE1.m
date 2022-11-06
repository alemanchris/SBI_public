%{
 EXAMPLE 1:
 Uses example data to run SBI
%}
clear 
close all
clc

path = cd;
% Load Data
load([path,'/input/exampledata'])
%% Required inputs
% Asign Region Time Series
idta(:,1) = data(:,3);
idta(:,2) = data(:,1);
% Set Time Variable
itm = data_time;
% Set policy implementation year/date
itp = 1960;
%% Optional inputs
% Region names
irnam = ['REST';'REG1']; % Write REG1 not REG_1, number of characters must be equal
% Name your Outcome variable
ionam = ['My Outcome'];
% Name the units your time variable ('year','day','time')
itnam = ['Year'];
% Mapping (2:linear or 3:quadratic)
inm = 2;
% Smoothing step 
ismo = [];
% Smoother
itsmo = [];
% Boostrap % 1:Boostrao(Defaulr) 0: Skip Boostrap Step
ib = []; 
%{
     0: Moving average, with windown = 9;
     1: Interpolation smoother CSAPS (DEFAULT) with smoothing parameter = 0.0008
     2: Chebyshev nodes with Cheby Regression
     3: B-Splines
     4: HP Filter (lambda = 50)
%}
SBI_nrm(idta,itm,itp,inm,ismo,itsmo,[],ib,irnam,itnam,ionam);








