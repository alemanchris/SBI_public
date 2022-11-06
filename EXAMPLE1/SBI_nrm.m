%{

Copyright (C) 2022 Christian Aleman, Christopher Busch, Alexander Ludwig, Raül Santaeulalia-Llopis
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation

Use this program at you own risk, please report any bugs to 
christian.c.aleman@gmail.com

DESCRIPTION
SBI_norm calls functions:
    -smooth_ts: Smooths series, using a chosen smoother
    -norm_func: Performs normalization C to T and T to C, computes policy estimates and
        figures, see, figure index at the end.
    -bb_perform: Boostraps results from norm_func
OPTPUT OF THE CODE:
    - Table of results with bootstrap 90% CI
    - Normalization figures


User can choose between linear nm = 2 or quadratic nm = 3 time stage
transform

INPUTS: 
    Required inputs
        -idta: Outcome Variable of two regions, Tx2 Matrix, where T (time observations e.g
        year) 
        -itm: Time units Tx1 (Years, Months Days, etc)
        -itp: Policy date (date after which one expects policy effects to be non-zero)
    Optional inputs
        Normalization options: 
        -inm:   2: (Default)Time stage transform 3:quadratic time stage transform 
        -ismo:  1: (Default)Smoothing 0: No smoothing
        -itsmo: Type of smoother
        -ici:   CI bands threshold (DEFAULT 0.1 that is 90% Bands)
            (you can change the default values manually, in section "Options for smoothing" )
            0: Moving average, with windown = 9;
            1: Interpolation smoother CSAPS (DEFAULT) with smoothing parameter = 0.0008
            2: Chebyshev nodes with Cheby Regression
            3: B-Splines
            4: HP filter (lambda=50)
        -ib : 0: No bootstrap 1: Boostrap (DEFAULT)
        Labels
        -irnam: Region labels
        -itnam: Time label
        -ionam: Outcome label

OUTPUT: 
    -psi_esti: Nx3 Mapping parameters 
    -pol_esti: Nx1 Policy Estimate in log difference
    -tp_norm1: Nx1 Location of policy date after normalization
    -fp_norm:  Nx1 Location of the first point in C to time of T
    -l_opt:    Nx1 Lenght of the overlap interval 
    -flag :    Nx1 1:Mapping didnt converge 2:Mapping delivers non-sence
-Figures Index:  
  Mapping C to T
    fig_1: Mapping function
    fig_2: Time Series Before Normalization
    fig_3: Time Series After Normalization
    fig_4: Zoom Overlap Interval Non Interpolated data
    fig_5: Zoom Overlap Interval + Interpolated data
    fig_6: Cum Gamma zoom Overlap Interval, Interpolated data
    fig_7: Cum Gamma Zoom Overlap Interval, + Non interpolated data
    fig_8: Combined figure 5 and figure 6
  Mapping T to C
    fig_9:  Time Series After Normalization
    fig_10: Zoom Overlap Interval Interpolated data
    fig_11: Cum Gamma Zoom Overlap Interval Interpolated data
    fig_12: Combined figure 10 and figure 11
  Confidence Bands
    fig_20: Cummulative Gamma, + median at point estimate
    fig_21: Time Series with bands
    fig_22: Time Series Median Boostrap draw with bands
    fig_23: histogram gamma around point estimate window
    fig_24: histogram gamma around median window
    fig_25: histogram gamma unrestricted
    fig_26: histogram gamma identification window length
    
    


Figure 1 has the option to have the mapping parameters in the
title, for this set opt.title_coeff=1

REMARKS: 
-Trimming(anchor) is optional, to trim select par.cut>0
-Local minizer is implemented

Before you run the this routine, make sure you have the following:
-At least two regions and their respective population
-Date of NATIONWIDE policy implementation. (itp)

Left to implement in this routine:
1. Input more than two regions


%}

function SBI_nrm(idta,itm,itp,inm,ismo,itsmo,ici,ib,irnam,itnam,ionam,varargin)
%clear
close all
clc

if nargin < 11
    ionam = [];
    if nargin < 10
        itnam = [];
        if nargin < 9
            irnam = [];
            if nargin < 8
                ib = [];
                if nargin < 7
                    ici = [];
                    if nargin < 6
                        itsmo = [];
                        if nargin < 5
                            ismo = [];
                            if nargin < 4
                                inm = [];
                                if nargin < 3
                                    error('User did not provide enough input arguments')
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Set paths and "names"
path = cd;
par.outpath = [path,'/output'];
par.inpath  = [path,'/input'];

%load([par.inpath,'/example_data']); 

% Set labels:
if isempty(ionam)
    par.outcome_name = ['outcome name'];
else
    par.outcome_name = ionam;
end
if isempty(itnam)
    par.time_name = ['time'];
else
    par.time_name = itnam;
end
if isempty(irnam)
    data.names_st = ['REG1';'REG2'];
else
    data.names_st = irnam; % ['REG01';'REG02';'REG03';'REG04';'REG05';'REG06';'REG07'];
end

% Set date of policy implementation: (must be in the same units as par.time)
%***************
par.tp = itp;
%***************

%% Data Treatment
% Extract Data
data.outcome = idta;
par.time = itm;       % This could be in Years, or Dates

par.ns = size(data.outcome,2);       % Number of regions 
if par.ns ==2
    par.nsp = 1;                  %
else
    % Input more than to regions, yet to be implemented
    error('More than 2 regions provided')
    par.nsp = par.ns+1;           % Number of regions + (1) Benchmark Region
end
par.nt = size(data.outcome,1);    % Number of time obeservations

%% Options for smoothing:
% Smoothing option
if isempty(ismo)
    opt.smooth = 1;                % 0: No smoothing 1: Smooth only-pre 2: Smoothing All(pre-post)
else
    opt.smooth = ismo;
end
% Time stage transform
if isempty(inm)
    par.nm = 2;           % 2: Linear mapping 3: Quadratic
else
    par.nm = inm;         % 2: Linear mapping 3: Quadratic
end
% For CSAPS
par.sp = 0.0008;%0.003;        % Smoothing parameter csaps
% For Moving average
opt.maw = 5;                   % Moving average window.
% For HP filer
opt.hpt = 50;                  % HP lambda
% For B-splines and Chebyshev Polinomial
opt.cheb_nodes = 1;            % (only for B-splines)1: use cheby nodes 2: Evenly Spaced 0: Automatic Matlab
par.m = 3;%round(0.9*par.nt);  % Mumber of collocation points: min: n+1,
par.n = 4;                     % Degree of the polinomial;

if isempty(itsmo)
    opt.type_sth = 1;  
else
    opt.type_sth = itsmo;  
end

%{
 0: Moving average
 1: Interpolation smoother CSAPS (DEFAULT)
 2: Chebyshev nodes with Cheby Regression
 3: B-Splines  
%}
%% Manual options for SBI Normalization:

opt.manual_guess = 2;   % 1: Pick Coefficient Initial Guess Manually 2: Polinomial analitical guess
par.m_guess = [0.8,2.14,0.75]; % Choose your manual guess, not recomended
if opt.manual_guess==1
    if size(par.m_guess,2)~=par.nm+1
        error('ERROR: Wrong number of normalization params provided for the guess')
    end
end
opt.wght = 1;           % Weighted estimation
opt.imeth = 'linear';   % Type of interpolation
par.cut = 0;            % Trimming, keep it at 0
par.cuts = 0;           % smooth starting in this cut
opt.sgues = 0;          % Use fmincon with sgues
par.min_match = 5;      % Minimum number of matched points
opt.graph_smooth = 1;   % Graph smoothed functions
opt.warn = 1;           % Show warnings for benchmark, warnings for boostrap automatically off
%% Boothstrapping params
if isempty(ib)
    opt.b = 1;         % 1: Bootstrap (automatically selects Block Bootstrapp)
else
    opt.b = ib; 
end
if isempty(ici)
    par.lev = 0.1; % 90%  
else
    par.lev = ici; % 90% 
end
par.nb = 3000;
par.nb_lim = 500;   % Stop when nb_lim number of boostrap iterations are sucessful
par.wdth = 4;       % length of bootstrap interval in Block Bootstrapping
par.iter = 0;       % initialize iteration counter, point estimate is iteration 0

%% Options for solvers 
opt.sol = optimset('Display','off','TolFun',1.0e-08,'TolX',1.0e-08);
opt.sol_fmin = optimoptions('fmincon','Display','off','FunctionTolerance',1.0e-10,'OptimalityTolerance',1.0e-10);
opt.optmin = optimoptions('fminunc','Display','off'); % for cheby


% Generate some parameters
par.t0  = par.time(1);
par.t00 = par.time(1); % Year where graph starts
par.t1 = par.time(end);

% Parameters for figures
par.wd = 0.2; % gamma restricted to wd% of the identification window
par.minyear = par.t0-1;     % Min year to graph
par.maxyear = par.t1;       % Max year to graph
%par.tu = par.tp-par.t0+1;   % Position of policy date
I = par.time<=par.tp;
par.tu = sum(I);   % Position of policy date
par.vis_ind = 0;   % 0: Turn off Benchmark Graphs (boostrap graphs still show) 1: Show all graphs
opt.title_coeff = 1;        % Add coefficients and flag to the title o fig1
opt.lw = 1.5;               % Linewidth
opt.lw2 = 1.7;
opt.fz2 = 19;% font size
opt.fz1 = 13;% font size
opt.fz3 = 23;% font size

% Initialize output vars

    olpCT = NaN(par.nsp,1);
    pol_estiCT = NaN(par.nsp,1);
    LS_CT = NaN(par.nsp,1);
    psi_estiCT = NaN(par.nsp,par.nm+1);
    tp_normCT  = NaN(par.nsp,1);
    fp_norm  = NaN(par.nsp,1);
    eflagCT  = NaN(par.nsp,1);

    olpTC = NaN(par.nsp,1);
    pol_estiTC = NaN(par.nsp,1);
    LS_TC = NaN(par.nsp,1);
    psi_estiTC = NaN(par.nsp,par.nm+1);
    tp_normTC  = NaN(par.nsp,1);
    eflagTC  = NaN(par.nsp,1);
    
    %% Determine Control and Treatment 
    flagC = detCT(data,par,opt);
    if flagC==1
        flagT = 2;
    else
        flagT = 1;
    end
    %% Define Control and Treatment
    %par.Cname = char(data.names_st(flagC,:));
    %par.Tname = char(data.names_st(flagT,:));
    
    par.Cname = ['$y_{\mathcal{C}}(t):$ ',char(data.names_st(flagC,:))];
    par.Tname = ['$y_{\mathcal{T}}(t):$ ',char(data.names_st(flagT,:))];

    par.Tnamenorm = ['$y_{\mathcal{T}}(t):$ ',char(data.names_st(flagT,:))];
    par.Cnamenorm = ['$y_{\mathcal{C}}(t;\mbox{\boldmath$\psi^{*}$}):$ ',char(data.names_st(flagC,:))];

    
    data.C = data.outcome(:,flagC);
    data.T = data.outcome(:,flagT);
    
    data.oriC = data.C;
    data.oriT = data.T;
    %% Smoothing
    [datas]=smooth_ts(data,par,opt);
    %% Perform SBI
      
    data.T = datas.T;
    data.C = datas.C;
    data.TO = datas.TO;
    data.CO = datas.CO;
    par.fapp = ['reg_',num2str(1),'sm_',num2str(opt.smooth),'typ_',num2str(opt.type_sth),'deg_',num2str(par.n),'knots_',num2str(par.m),'bb_',num2str(opt.b),'wght_',num2str(opt.wght),'quad_',num2str(par.nm)];
    par.name_fapp = 'nb500'; % for temporary resutls and boostrap figs;

    
        % Graph Smooth
    if opt.graph_smooth==1
        figure(1000)
        hold on
        plot(par.time,datas.T,'r-')
        plot(par.time,datas.C,'b-')
        plot(par.time,datas.TO,'rx')
        plot(par.time,datas.CO,'bo')
    end
    
    [out_r] = norm_func(data,par,opt);
    
    
    %% Extract  Results
    i = 1; % if i>1, means there are more regions, yet to implement
    olpCT(i,1) = out_r.olpCT;
    olpTC(i,1) = out_r.olpTC;
    pol_estiCT(i,1) = out_r.pol_estiCT;
    pol_estiTC(i,1) = out_r.pol_estiTC;    
    psi_estiCT(i,:) = out_r.psi_estiCT;       
    psi_estiTC(i,:) = out_r.psi_estiTC;
    LS_CT(i,1) = out_r.LS_CT;       
    LS_TC(i,1) = out_r.LS_TC;
    tp_normCT(i,1)  = out_r.tp_normCT;
    tp_normTC(i,1)  = out_r.tp_normTC;
    fp_norm(i,1)  = out_r.fp_norm;

    
    %% Perform bootstrap
    
    bb_perform(datas,data,out_r,par,opt);

    
    
  %% Save Results:
    % Currently the routine saves the results from the display table, this
    % is done within bb_perform
    %{
    fname=['output/olpCT_', par.fapp, '.txt'];
    save(fname,'olpCT','-ascii','-double','-tabs');

    fname=['output/psi_estiCT_', par.fapp, '.txt'];
    save(fname,'psi_estiCT','-ascii','-double','-tabs');
    
    fname=['output/LS_CT_', par.fapp, '.txt'];
    save(fname,'LS_CT','-ascii','-double','-tabs');
    
    fname=['output/pol_estiCT_', par.fapp, '.txt'];
    save(fname,'pol_estiCT','-ascii','-double','-tabs');
    
    fname=['output/fp_norm', par.fapp, '.txt'];
    save(fname,'fp_norm','-ascii','-double','-tabs');
    
    fname=['output/tp_normCT_', par.fapp, '.txt'];
    save(fname,'tp_normCT','-ascii','-double','-tabs');
    
    fname=['output/olpTC_', par.fapp, '.txt'];
    save(fname,'olpTC','-ascii','-double','-tabs');

    fname=['output/psi_estiTC_', par.fapp, '.txt'];
    save(fname,'psi_estiTC','-ascii','-double','-tabs');
    
    fname=['output/LS_TC_', par.fapp, '.txt'];
    save(fname,'LS_TC','-ascii','-double','-tabs');
    
    fname=['output/pol_estiTC_', par.fapp, '.txt'];
    save(fname,'pol_estiTC','-ascii','-double','-tabs');
    
    fname=['output/tp_normTC_', par.fapp, '.txt'];
    save(fname,'tp_normTC','-ascii','-double','-tabs');
    %}   


disp('End of Routine')
end
%--------------------------------------------------------------------------
function [datas,ind_smo]=smooth_ts(data,par,opt)
%{
Input:
    data.outcome: Raw Series
Output:
    data.soutcome: Smooth Outcome

Type of smoother
 0: Moving average
 1: Interpolation smoother CSAPS (DEFAULT)
 2: Chebyshev nodes with Cheby Regression
 3: B-Splines  
%}
% Keep original data:
    ind_smo = 0; % 1: Smoothing step failed (deprecated)
    datas.TO = data.T;
    datas.CO = data.C;
    datas.T = data.T;
    datas.C = data.C;

    if opt.smooth==1
        pos = par.tu;
    elseif opt.smooth==2
        pos = par.nt;
    end

    if opt.smooth~=0
       % Smooth only positive values       
       Taux = data.T(1:pos);
       IT = Taux>par.cuts;       
       [~,AA] = max(IT);
       IT(AA:end) = true;
       nT = sum(IT);
       Caux = data.C(1:pos);
       IC = Caux>par.cuts;
       [~,AA] = max(IC);
       IC(AA:end) = true;
       nC = sum(IC);
       try
           if opt.type_sth==1     % Poynomial Smoother (DEFAULT)
               datas.T(IT) = csaps((1:nT)',Taux(IT),par.sp,(1:nT)');
               datas.C(IC) = csaps((1:nC)',Caux(IC),par.sp,(1:nC)');

           elseif opt.type_sth==2 %  Cheby Polynomials with Cheby nodes
               datas.T(IT) = coll((1:nT)',Taux(IT),par,opt);
               datas.C(IC) = coll((1:nC)',Caux(IC),par,opt);
           elseif opt.type_sth==3 % B-Splines 
               datas.T(IT) = spline_reg((1:nT)',Taux(IT),par,opt);
               datas.C(IC) = spline_reg((1:nC)',Caux(IC),par,opt);
           elseif opt.type_sth==0   % Moving average 
               datas.T(1:pos) = smoothdata(data.T(1:pos),'movmean',opt.maw);
               datas.C(1:pos) = smoothdata(data.C(1:pos),'movmean',opt.maw);
           elseif opt.type_sth==4   % HP filter 
               datas.T(1:pos) = hpfilter(data.T(1:pos),opt.hpt); 
               datas.C(1:pos) = hpfilter(data.C(1:pos),opt.hpt); 
              
           end
       catch
           error('Error While Smoothing');
       end
    end
end
%-------------------------------------------------------------------------
function [o_vals]=norm_func(data,par,opt)
close all


%{
figure(100)
hold on
plot(data.CO,'bx')
plot(data.C,'b--')
plot(data.TO,'rx')
plot(data.T,'r--')
%}



% Keep full series
par.nt_long = length(data.C);
par.time_long = par.time;

% Trim to use prepolicy data for normalization
par.time = par.time(1:par.tu,1);
data.C = data.C(1:par.tu,1);

par.nt = length(data.C);

% Recompute some location parameters according to the trimming
data.T = data.T(1:par.tu);


par.tvec_c      = (1:par.nt)';
par.tvec_c_long = (1:par.nt_long)';



% Initial guess for mapping parameters
%{ 
phi(1): Magnitude Shifter
phi(2): Time Shifter
phi(3): Speed Shifter
%}
%% Mapping C to T
m_guess = par.m_guess;

% Activar cuando es a manubrio
%{
opt.CT = 1;
[cf.phiCT,x0] = make_x(data,opt,m_guess,par); 
func_df(x0,data.T,data.C,opt,par)
dert_debug = 1;
%}
opt.CT = 1;
[cf.phiCT,x0] = make_x(data,opt,m_guess,par); 
if opt.manual_guess==1
    if size(cf.phiCT,2)~=par.nm+1
        error('ERROR: The number of manual normalization params, is not the correct one')
    end
end
if opt.sgues==1
    [x0,~,~] = fmincon(@(x)func_df(x,data.T,data.C,opt,par),x0,[],[],[],[],[0,-Inf,-Inf],[Inf,Inf,Inf],[],opt.sol_fmin);
end
try
    [xCT,fvalCT,eflagCT] = fminunc(@func_df,x0,opt.sol,data.T,data.C,opt,par);
catch
    xCT = m_guess;
    fvalCT = 0;
    eflagCT = -3; 
end
if eflagCT==0
    if opt.warn ==1
    disp('WARNING! SOLVER C to T QUIT DUE TO MAX ITER')
    end
end
cf.phiCT = xCT;


[~,~,nmatched]  = func_df(xCT,data.T,data.C,opt,par);
if nmatched<par.min_match
    if opt.warn ==1
    disp('Not Enough matched points')    
    end
end


tu_norm = func_stage(par.tu,cf.phiCT,0,par);
tp_normCT_1 = interp1(par.tvec_c_long,par.time_long,tu_norm); % Convert to time units
aux1=func_stage(par.tu,cf.phiCT,1,par);
tp_normCT_2 = interp1(par.tvec_c_long,par.time_long,aux1);
%% Mapping T to C
opt.manual_guess = 1; 
if par.nm == 3
    m_guess = [1/xCT(1),-xCT(2)/xCT(3) ,1/xCT(3),0];
else
    m_guess = [1/xCT(1),-xCT(2)/xCT(3) ,1/xCT(3)];
end
opt.CT = 0;
[cf.phiTC,x0] = make_x(data,opt,m_guess,par); 
if opt.sgues ==1
   [x0,~,~] = fmincon(@(x)func_df(x,data.C,data.T,opt,par),x0,[],[],[],[],[0,-Inf,-Inf],[Inf,Inf,Inf],[],opt.sol_fmin);
end
try 
    [xTC,fvalTC,eflagTC] = fminunc(@func_df,x0,opt.sol,data.C,data.T,opt,par);
catch
    xTC = x0;
    fvalTC = 0;
    eflagTC = -3; 
end
if eflagTC==0
    if opt.warn ==1
    disp('WARNING! SOLVER T to C QUIT DUE TO MAX ITER')
    end
end
cf.phiTC = xTC;


aux1=func_stage(par.tu,cf.phiTC,0,par);
tp_normTC_1 = interp1(par.tvec_c_long,par.time_long,aux1);
aux1=func_stage(par.tu,cf.phiTC,1,par);
tp_normTC_2 = interp1(par.tvec_c_long,par.time_long,aux1);
aux1=func_stage(par.nt_long,cf.phiTC,0,par);
tp_normTC_A = interp1(par.tvec_c_long,par.time_long,aux1);

aux1=func_stage(1,cf.phiCT,0,par);
yr_init = interp1(par.tvec_c_long,par.time_long,aux1,'linear',NaN);
if par.iter==0
disp(['Coefficients in numerical mapping C to T: ', num2str(cf.phiCT)]);
disp(['Coefficients in numerical mapping T to C: ', num2str(cf.phiTC)]);
disp('Results from C to T:')
disp(['Year in ',par.Cname,': ', num2str(par.tp), ' Stage in ',par.Tname,': ', num2str(tp_normCT_1)]);
disp(['Year in ',par.Cname,': ', num2str(par.t0), ' Stage in ',par.Tname,': ', num2str(yr_init)]);
disp('Results from T to C:')
disp(['Year in ',par.Tname,': ', num2str(par.tp), ' Stage in ',par.Cname,': ', num2str(tp_normTC_1)]);
disp(['Year in ',par.Tname,': ', num2str(par.t1), ' Stage in ',par.Cname,': ', num2str(tp_normTC_A)]);
end
%% Compute Relevant Output:
% Normalized series
svec = par.tvec_c_long;
tvec = func_stage(svec,cf.phiCT(:),0,par);    % corresponding time in ROECD (i.e., when will ROECD look like Germany)

time_long_y = interp1(par.tvec_c_long,par.time_long,tvec,'linear',NaN); % Convert to time units
C_norm = cf.phiCT(1).*interp1(tvec(1:par.tu),data.C,par.tvec_c_long,'linear',NaN);
CO_norm = cf.phiCT(1).*interp1(tvec,data.CO,par.tvec_c_long,'linear',NaN);
C_norm_no = cf.phiCT(1).*data.CO;
T_interp = interp1(par.tvec_c_long,data.TO,tvec,'linear',NaN);

tvec_TC = func_stage(svec,cf.phiTC(:),0,par);    % corresponding time in ROECD (i.e., when will ROECD look like Germany)
T_norm_TC = cf.phiTC(1).*interp1(tvec_TC(1:par.tu),data.T,par.tvec_c_long,'linear',NaN);
TO_norm_TC = cf.phiTC(1).*interp1(tvec_TC,data.TO,par.tvec_c_long,'linear',NaN);



%{
figure(1)
hold on
plot(par.tvec_c_long,data.CO,'b-')
plot(par.tvec_c_long,CO_norm,'b--')
plot(par.tvec_c_long,data.TO,'r-')
%
figure(2)
hold on
plot(par.tvec_c_long,data.CO,'b-')
plot(par.tvec_c_long,TO_norm,'r--')
plot(par.tvec_c_long,data.TO,'r-')
%
figure(3)
hold on
plot(par.tvec_c_long,data.CO,'b-')
plot(par.tvec_c_long,CO_norm_TC,'b--')
plot(par.tvec_c_long,data.TO,'r-')
%
figure(4)
hold on
plot(par.tvec_c_long,data.CO,'b-')
plot(par.tvec_c_long,TO_norm_TC,'r--')
plot(par.tvec_c_long,data.TO,'r-')
%}

%% Policy Effect CT
In = sum(isnan(CO_norm([par.tu-1,par.tu,par.tu+1])));

pdiff_CT = NaN(size(data.TO));
pdiff_CT_int = NaN(size(data.TO));
pdiff_CT_pre = NaN(size(data.TO)); % Prepolicy gamma
pdiff_CT_int_pre = NaN(size(data.TO));% Prepolicy interpolated gamma
tp_normCT_3=func_stage(par.tu,cf.phiCT,0,par);
if tp_normCT_3>size(data.TO,1)
    % Extreme case
    if opt.warn ==1
    disp('ERROR: Detected extreme case CT')
    end
    tp_normCT_3 = size(data.TO,1)-1;
    In = 1;
end

% Flag to fish for outliers
flag_out = 1; % 1: outlier(flip etc...)
flag_out_flip = 0;
flag_out_sign = 0;
flag_out_nons = 0;
if In==0 && abs(cf.phiCT(2))<5000

  %% Estimated Policy Effect
     

     % For interpolated T
     Iaux = (tvec - par.tu)>0; % not larger or equal cuz if equal I the effect of the first day =0
     [~,fpos] = max(Iaux); % Position of first day
 
     dataX = [tvec(fpos:par.tu);par.tvec_c_long(par.tu:floor(tp_normCT_3))];
     dataYr = [T_interp(fpos:par.tu);data.TO(par.tu:floor(tp_normCT_3))];
     dataYb = [C_norm_no(fpos:par.tu);CO_norm(par.tu:floor(tp_normCT_3))];
     aux_data=sortrows([dataX,dataYr,dataYb]);
     size_aux = size(dataX,1);
     pdiff_CT_alt = NaN(size_aux,1); 
     for ii = 2:size_aux
         unr = trapz(aux_data(1:ii,1)',aux_data(1:ii,2)'); % Area under red
         unb = trapz(aux_data(1:ii,1),aux_data(1:ii,3));   % Area under blue
         pdiff_CT_alt(ii) = (unr-unb)/unb;
         % Lets asign it the position it deserves for the graph
         I = aux_data(ii,1)==tvec;
         if sum(I)==1
            pdiff_CT_int(I)=pdiff_CT_alt(ii);
         else
            I =  aux_data(ii,1)==par.tvec_c_long;
            pdiff_CT(I)=pdiff_CT_alt(ii);
         end
     end
     olp_CT = tp_normCT_1-par.tp;
     if olp_CT<=0
         gammarCT = NaN;
         LS_CT = NaN;
         DLS_CT = NaN; % Stagely
         ind_del = 1;
     else
         gammarCT = pdiff_CT_alt(end);
         LS_CT = unr-unb;
         DLS_CT = LS_CT/olp_CT; % Stagely
      

         %% Pre-policy (Not Gamma)

         pinit_aux = (par.tu-floor(olp_CT)); 
         fpos_1 = (fpos-1);
         I_aux = (tvec-pinit_aux)>0;
         [~,fpos_2] = max(I_aux);

         dataX = [tvec(fpos_2:fpos_1);par.tvec_c_long(pinit_aux:par.tu)];
         dataYr = [T_interp(fpos_2:fpos_1);data.TO(pinit_aux:par.tu)];
         dataYb = [C_norm_no(fpos_2:fpos_1);CO_norm(pinit_aux:par.tu)];
         aux_data=sortrows([dataX,dataYr,dataYb]);
         size_aux = size(dataX,1);
         pdiff_CT_alt = NaN(size_aux,1); 
         for ii = 2:size_aux
             unr = trapz(aux_data(1:ii,1)',aux_data(1:ii,2)'); % Area under red
             unb = trapz(aux_data(1:ii,1),aux_data(1:ii,3)); % Area under blue
             pdiff_CT_alt(ii) = (unr-unb)/unb;
             % Lets asign it the position it deserves for the graph
             I = aux_data(ii,1)==tvec;
             if sum(I)==1
                pdiff_CT_int_pre(I)=pdiff_CT_alt(ii);
             else
                I =  aux_data(ii,1)==par.tvec_c_long;
                pdiff_CT_pre(I)=pdiff_CT_alt(ii);
             end
         end
     end
      % Normalize fvalCT
     
     aux_Cnorm = cf.phiCT(1).*interp1(tvec(1:par.tu),data.C(1:par.tu),svec(1:par.tu),'linear',NaN);     
     avT = mean(data.TO(1:par.tu)); % Average deaths in reference region
     I = not(isnan(aux_Cnorm));   
     numpoCT = sum(I);  % Number of matcheed points
     fvalCT =  fvalCT./(avT.*numpoCT); % Normalize it
     
     flag_out = 0;
     %{
     Nonsence can be:
     -a flip
     -extreme cases
     -non-convergence
     -crazy psis
     -gamma = NaN
     %}
     % Fish for outliers
     if isnan(gammarCT)
         if ind_del ~= 1
            flag_out = 1; % 
            flag_out_nons = 1; % Flag nonsence
         end
     end
     if olp_CT<0
        flag_out = 1;  % Flip
        flag_out_flip = 1;
     end      
     if LS_CT>=0
        %flag_out = 1; % I should not restric these.
        flag_out_sign = 1;
     end
     if eflagCT<=0
       flag_out = 1; % Non Convergence, or corner solution
       flag_out_nons = 1;
     end  
else
    %I = not(isnan(ln_gdpn_rest_sc));

     gammarCT = NaN;
     olp_CT = NaN;
     LS_CT = NaN;
     DLS_CT = NaN;    
    
end

%% Policy Effect TC
%In = sum(isnan(TO_norm_TC([par.tu-1,par.tu,par.tu+1])));
In = sum(isnan(TO_norm_TC([par.tu-1,par.tu]))); % Otherwise I cannot compute policy effect
% Same as gamma
pdiff_TC = NaN(size(data.TO));
tp_normTC_3=func_stage(par.tu,cf.phiTC,0,par);
if tp_normTC_3<1
    % Extreme case
    if opt.warn ==1
    disp('ERROR: Detected extreme case TC')
    end
    tp_normTC_3 = 2;
end
if In==0 && abs(cf.phiTC(2))<5000

    %% Estimated Policy Effect
  
     auxTO_norm = cf.phiTC(1).*data.TO;
     Cpol = interp1(par.tvec_c_long,data.CO,tp_normTC_3,'linear');
     ldr =  auxTO_norm(par.tu)*(ceil(tp_normTC_3)-tp_normTC_3);% Last day red
     ldb =  Cpol*(ceil(tp_normTC_3)-tp_normTC_3);% Last day blue

     for ij = ceil(tp_normTC_3):par.tu
        unr=sum(TO_norm_TC(ceil(tp_normTC_3):ij))+ldr; % Area under Blue
        unb=sum(data.CO(ceil(tp_normTC_3):ij))+ldb;
        pdiff_TC(ij) = (unr-unb)/unb;
     end
     unr=sum(TO_norm_TC(ceil(tp_normTC_3):par.tu))+ldr; % Area under Blue
     unb=sum(data.CO(ceil(tp_normTC_3):par.tu))+ldb;
     gammarTC = (unr-unb)/unb;
     olp_TC = par.tp-tp_normTC_1;
     LS_TC = unr-unb;
     DLS_TC = LS_TC/olp_TC;
    
     aux_Tnorm = cf.phiTC(1).*interp1(tvec_TC(1:par.tu),data.T(1:par.tu),svec(1:par.tu),'linear',NaN);     
     avC = mean(data.CO(1:par.tu)); % Average deaths in reference region
     I = not(isnan(aux_Tnorm));   
     numpoTC = sum(I);  % Number of matched points    
     fvalTC =  fvalTC./(avC.*numpoTC); % Normalize it:
     
else
    %I = not(isnan(ln_gdpn_rest_sc));

     gammarTC = NaN;
     olp_TC = NaN;
     LS_TC = NaN;
     DLS_TC = NaN;    
    
end

% Display Policy effects:
if par.iter==0
disp(['Policy Effect C to T: ',num2str(round(gammarCT,3))])
disp(['Policy Effect T to C: ',num2str(round(gammarTC,3))])
disp(['fvalCT  : ',num2str(fvalCT)])
disp(['fvalTC  : ',num2str(fvalTC)])
end



II = eflagCT==-1||eflagCT==-3; % Algo stalled 
III = isnan(tp_normCT_1); % Extreme case
if III==1
    o_vals.flagCT = 'Extreme Case';
end
if II ==1
    o_vals.flagCT = 'Algorithm Didnt Converge';
end

if  II==1 || III==1
    if opt.warn ==1
    disp(['Mapping CT Not Working for region: ',par.Cname,' Flag: ',o_vals.flagCT])
    end
    o_vals.olpCT = -999;
    o_vals.pol_estiCT = -999;
    o_vals.LS_CT = -999;  
    o_vals.psi_estiCT = [-999,-999,-999];
    o_vals.tp_normCT  = -999;
    o_vals.fp_norm  = -999;
    o_vals.fvalCT =  -999;
    o_vals.flag_out = 1;
    o_vals.flag_out_nons = 1;
    o_vals.flag_out_flip = 0;
    o_vals.flag_out_sign = 0;
    o_vals.pdiff_CT = NaN(size(pdiff_CT));
    o_vals.pdiff_CT_int= NaN(size(pdiff_CT));
    o_vals.pdiff_CT_pre= NaN(size(pdiff_CT));
    o_vals.pdiff_CT_int_pre= NaN(size(pdiff_CT));  

    
    o_vals.TO = ones(size(data.CO)).*(-999);
    o_vals.CO = ones(size(data.CO)).*(-999);  
    o_vals.T = ones(size(data.CO)).*(-999);
    o_vals.C = ones(size(data.CO)).*(-999);  
    o_vals.Cnorm = ones(size(data.CO)).*(-999);  
    o_vals.CO_norm = ones(size(data.CO)).*(-999);
    o_vals.time = par.time_long;
    o_vals.time_norm = par.time_long;
    
    o_vals.fvalCT =  -999;
    o_vals.flag_out = 1; % 
    o_vals.flag_out_nons = 1;
    o_vals.flag_out_flip = 0;
    o_vals.flag_out_sign = 0;
    
else     

    o_vals.olpCT = olp_CT;
    o_vals.pol_estiCT = gammarCT;  
    o_vals.LS_CT = LS_CT;  
    o_vals.psi_estiCT = cf.phiCT;
    o_vals.tp_normCT  = tp_normCT_1;
    o_vals.fp_norm  = yr_init;
    o_vals.fvalCT =  fvalCT;
    o_vals.flag_out = flag_out;
    o_vals.flag_out_nons = flag_out_nons;
    o_vals.flag_out_flip = flag_out_flip;
    o_vals.flag_out_sign = flag_out_sign;
    
    o_vals.pdiff_CT = pdiff_CT;
    o_vals.pdiff_CT_int = pdiff_CT_int;
    o_vals.pdiff_CT_pre = pdiff_CT_pre;
    o_vals.pdiff_CT_int_pre = pdiff_CT_int_pre;
    
    o_vals.TO = data.TO;
    o_vals.CO = data.CO;  
    o_vals.CO_norm = CO_norm;
    o_vals.T = data.T;
    o_vals.C = data.C;  
    o_vals.Cnorm = C_norm;  
    o_vals.time = par.time_long;
    o_vals.time_norm = time_long_y;
end


    
II = eflagTC==-1||eflagTC==-3;
III = isnan(tp_normTC_1);

if  III==1
    o_vals.flagTC = 'Extreme Case';
end
if II ==1
    o_vals.flagTC = 'Algorithm Didnt Converge';
end

if  II==1 || III==1
    if opt.warn ==1
    disp(['Mapping TC Not Working for region: ',par.Cname,' Flag: ',o_vals.flagTC])
    end
    o_vals.olpTC = -999;
    o_vals.pol_estiTC = -999;
    o_vals.LS_TC = -999;  
    o_vals.psi_estiTC = [-999,-999,-999];
    o_vals.tp_normTC  = -999;
    o_vals.pdiff_TC = ones(size(pdiff_TC)).*(-999);
    
else     


    o_vals.olpTC = olp_TC;   %par.tp-tp_normTC_1;
    o_vals.pol_estiTC = gammarTC;
    o_vals.LS_TC = LS_TC;  
    o_vals.psi_estiTC = cf.phiTC;
    o_vals.tp_normTC  = tp_normTC_1;
    o_vals.fvalTC =  fvalTC; 
    o_vals.pdiff_TC = pdiff_TC;

end



%% Generate and save Figures:
%
if par.vis_ind ==1
par.visible = 'on';        % 'on' if you want to see the graphs as they come
if isnan(tp_normTC_1)
    tp_normTC_1 = par.t0;
    tp_normTC_3 = par.t0;
end

if isnan(tp_normCT_1)
    tp_normCT_1 = par.t1;
    tp_normCT_3 = par.t1;
end



% Generate title:
aosCT = sprintf('%.2f,' , [cf.phiCT]);
aosCT = aosCT(1:end-1);  % strip final comma

aosTC = sprintf('%.2f,' , [cf.phiTC]);
aosTC = aosTC(1:end-1);  % strip final comma

%% Mapping Function
f1=figure(1);
set(f1, 'Visible',par.visible);
hold on
plot(svec,tvec,'k-','linewidth',1.1)
plot(svec,svec,'k--','linewidth',1)
legend(par.Cname,'Location','Best','interpreter','latex')
xlabel('Time ')
ylabel('Stage')
title(aosCT);
print('-depsc',['output/f1_', par.fapp, '.eps']);
saveas(f1,['output/f1_', par.fapp,'.png']);
hold off

%% Before Normalization
ymin = min([min(data.TO),min(data.CO),min(CO_norm)]);
ymax = max([max(data.TO),max(data.CO),max(CO_norm)])*1.1;
ymax2 = max([max(data.TO),max(CO_norm)])*1.1;


f2=figure(2);
set(f2, 'Visible',par.visible);
hold on
l1 = plot(par.time_long(1:end),data.TO,'r^','linewidth',opt.lw);
l2 = plot(par.time_long,data.CO,'bo','linewidth',opt.lw);
plot(par.time_long(1:par.tu),data.T(1:par.tu),'r-','linewidth',opt.lw)
plot(par.time_long(1:par.tu),data.C(1:par.tu),'b-','linewidth',opt.lw)
xline(par.tp,'k-','linewidth',1.1)
legend([l2,l1],{par.Cname,par.Tname},'Location','Best','interpreter','latex','FontSize',opt.fz2-1,'box','off')
xlim([par.t00,par.t1]);
ylim([ymin,ymax]);
set(gca, 'FontSize',opt.fz1);
ylabel(par.outcome_name ,'interpreter','latex','FontSize',opt.fz2)
xlabel(par.time_name,'interpreter','latex','FontSize',opt.fz2)
print('-depsc',['output/f2_', par.fapp, '.eps']);
saveas(f2,['output/f2_', par.fapp,'.png']);
hold off

%% After Normalization C to T 


f3=figure(3);
set(f3, 'Visible',par.visible);
hold on
%l2 = plot(time_long_y,T_interp,'r^','linewidth',opt.lw);
l2 = plot(par.time_long(1:end),data.TO,'r^','linewidth',opt.lw);
%l1 = plot(par.time_long,CO_norm,'bo','linewidth',opt.lw);
l1 = plot(time_long_y,C_norm_no,'bo','linewidth',opt.lw);
plot(par.time_long(1:par.tu),C_norm(1:par.tu),'b--','linewidth',opt.lw)
plot(par.time_long(1:par.tu),data.T(1:par.tu),'r-','linewidth',opt.lw)
area([par.tp,tp_normCT_1],[9999 9999],-9999,'FaceColor','r','FaceAlpha',.1,'EdgeColor','None','linestyle','none');
xline(par.tp,'k-','linewidth',1.1)
xline(tp_normCT_1,'k--','linewidth',1.1)
legend([l1,l2],{[par.Cnamenorm,' Norm'],par.Tnamenorm},'Location','Best','interpreter','latex','FontSize',opt.fz2-1,'box','off')
xlim([par.t00,par.t1]);
ylim([ymin,ymax2]);
set(gca, 'FontSize',opt.fz1);
ylabel(par.outcome_name,'interpreter','latex','FontSize',opt.fz2)
xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
%title(['Policy Effect: ',num2str(round(gammarCT,3)),' MinError: ',num2str(fvalCT)],'FontSize',opt.fz1-2)
%subtitle(['Phis: ',num2str(cf.phiCT(1)),', ',num2str(cf.phiCT(2)),', ',num2str(cf.phiCT(3))],'FontSize',opt.fz1-2)
print('-depsc',['output/f3_', par.fapp, '.eps']);
saveas(f3,['output/f3_', par.fapp, '.png']);
hold off

% Zoom 
if par.tp-tp_normCT_1>0
    % A flip
    error('Flip: You should not be here, check why!')

else
    if ceil(tp_normCT_3)>=par.nt_long
        tp_normCT_3 = ceil(tp_normCT_3)-1;
    end
    ymin = min([min(data.TO(par.tu-3:ceil(tp_normCT_3)+1)),min(CO_norm(par.tu-3:ceil(tp_normCT_3)+1))])*0.95;
    ymax = max([max(data.TO(par.tu-3:ceil(tp_normCT_3)+1)),max(CO_norm(par.tu-3:ceil(tp_normCT_3)+1))])*1.05;

    f4=figure(4);
    set(f4, 'Visible',par.visible);
    hold on
    l1 = plot(time_long_y(1:par.tu),T_interp(1:par.tu),'r^','MarkerSize',8,'linewidth',opt.lw2);
    %l1 = plot(par.time_long(1:end),data.TO,'r^','MarkerSize',8,'linewidth',opt.lw2);
    %l2 = plot(par.time_long,CO_norm,'bo','MarkerSize',8,'linewidth',opt.lw2);
    l2 = plot(time_long_y(1:par.tu),C_norm_no(1:par.tu),'bo','MarkerSize',8,'linewidth',opt.lw);
    area([par.tp,tp_normCT_1],[9999 9999],-9999,'FaceColor','r','FaceAlpha',.1,'EdgeColor','None','linestyle','none');
    xline(par.tp,'k-','linewidth',1.1)
    xline(tp_normCT_1,'k--','linewidth',1.1)
    legend([l2,l1], {[par.Cnamenorm,' Norm'],par.Tnamenorm},'Location','Best','interpreter','latex','FontSize',opt.fz2-1,'box','off')
    xlim([par.tp-3,tp_normCT_1+1]);
    ylim([ymin,ymax])
    set(gca, 'FontSize',opt.fz1);
    ylabel(par.outcome_name ,'interpreter','latex','FontSize',opt.fz2)
    xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
    print('-depsc',['output/f4_', par.fapp, '.eps']);
    saveas(f4,['output/f4_', par.fapp, '.png']);
    hold off

    f5=figure(5);
    set(f5, 'Visible',par.visible);
    hold on
    l0 = plot(time_long_y(1:par.tu),T_interp(1:par.tu),'ro','MarkerSize',8,'linewidth',opt.lw2);
    l1 = plot(par.time_long(1:floor(tu_norm)),data.TO(1:floor(tu_norm)),'r^','MarkerSize',8,'linewidth',opt.lw2);
    l3 = plot(par.time_long(1:floor(tu_norm)),CO_norm(1:floor(tu_norm)),'b^','MarkerSize',8,'linewidth',opt.lw2);
    l2 = plot(time_long_y(1:par.tu),C_norm_no(1:par.tu),'bo','MarkerSize',8,'linewidth',opt.lw);
    area([par.tp,tp_normCT_1],[9999 9999],-9999,'FaceColor','r','FaceAlpha',.1,'EdgeColor','None','linestyle','none');
    xline(par.tp,'k-','linewidth',1.1)
    xline(tp_normCT_1,'k--','linewidth',1.1)
    legend([l2,l1], {[par.Cnamenorm,' Norm'],par.Tnamenorm},'Location','Best','interpreter','latex','FontSize',opt.fz2-1,'box','off')
    xlim([par.tp-3,tp_normCT_1+1]);
    ylim([ymin,ymax])
    set(gca, 'FontSize',opt.fz1);
    ylabel(par.outcome_name ,'interpreter','latex','FontSize',opt.fz2)
    xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
    print('-depsc',['output/f5_', par.fapp, '.eps']);
    saveas(f5,['output/f5_', par.fapp, '.png']);
    hold off

    %% Policy Effect Zoom
    alt_yr = [par.tp,tp_normCT_1];
    sizealt=size(alt_yr,2);

    if o_vals.pol_estiCT<0
    ymin2 = min([min(pdiff_CT),o_vals.pol_estiCT])*1.05;
    ymax2 = max([max(pdiff_CT_pre),+eps,max(pdiff_CT)])*1.05;
        if isnan(ymin2)  
            ymin2 = -0.1; 
        end
    else
    ymin2 = min([min(pdiff_CT_pre),-eps,min(pdiff_CT)])*1.05;%-eps;
    ymax2 = max([max(pdiff_CT),o_vals.pol_estiCT])*1.05;
        if isnan(ymax2)  
            ymax2 = 0.1; 
        end
    end


    f6 = figure(6);
    set(f6, 'Visible',par.visible);
    hold on
    plot(alt_yr,repmat(o_vals.pol_estiCT,[1,sizealt]),'m-','linewidth',opt.lw+1)
    plot(par.time_long,pdiff_CT,'ko','MarkerFaceColor','k','MarkerSize',10)
    plot(tp_normCT_1,gammarCT,'ko','MarkerFaceColor','k','MarkerSize',10)
    plot(par.time_long,zeros(par.nt_long,1),'Color', [0 0 0]+0.2)
    area([par.tp,tp_normCT_1],[9999 9999],-9999,'FaceColor','r','FaceAlpha',.1,'EdgeColor','None','linestyle','none');
    xline(par.tp,'k-','linewidth',1.1)
    xline(tp_normCT_1,'k--','linewidth',1.1)
    %xlim([par.minyear,par.maxyear])
    xlim([par.tp-3,tp_normCT_1+1]);
    ylim([ymin2,ymax2])
    legend('$\gamma$','Location','northwest','FontSize',opt.fz2-1,'interpreter','latex','box','off')
    set(gca, 'FontSize',opt.fz1);
    ylabel('$\gamma(s)$','interpreter','latex','FontSize',opt.fz2)
    xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
    print('-depsc',['output/f6_', par.fapp, '.eps']);
    saveas(f6,['output/f6_', par.fapp, '.png']);
    hold off

    f7 = figure(7);
    set(f7, 'Visible',par.visible);
    hold on
    plot(alt_yr,repmat(o_vals.pol_estiCT,[1,sizealt]),'m-','linewidth',opt.lw+1)
    plot(par.time_long,pdiff_CT,'^','MarkerSize',10,'linewidth',2.2,'Color',[0.4940 0.1840 0.5560])
    plot(time_long_y,pdiff_CT_int,'o','MarkerSize',10,'linewidth',2.2,'Color',[0.4940 0.1840 0.5560])
    plot(par.time_long,pdiff_CT_pre,'k^','MarkerSize',10,'linewidth',2.2)
    plot(time_long_y,pdiff_CT_int_pre,'ko','MarkerSize',10,'linewidth',2.2)
    plot(tp_normCT_1,gammarCT,'ko','MarkerFaceColor','k','MarkerSize',10)
    plot(par.time_long,zeros(par.nt_long,1),'Color', [0 0 0]+0.2)
    area([par.tp,tp_normCT_1],[9999 9999],-9999,'FaceColor','r','FaceAlpha',.1,'EdgeColor','None','linestyle','none');
    xline(par.tp,'k-','linewidth',1.1)
    xline(tp_normCT_1,'k--','linewidth',1.1)
    %xlim([par.minyear,par.maxyear])
    xlim([par.tp-3,tp_normCT_1+1]);
    ylim([ymin2,ymax2])
    legend('$\gamma$','Location','northwest','FontSize',opt.fz2-1,'interpreter','latex','box','off')
    set(gca, 'FontSize',opt.fz1);
    ylabel('$\gamma(s)$','interpreter','latex','FontSize',opt.fz2)
    xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
    print('-depsc',['output/f7_', par.fapp, '.eps']);
    saveas(f7,['output/f7_', par.fapp, '.png']);
    hold off

    f8 = figure(8);
    set(f8, 'Visible',par.visible);
    subplot(2,1,1)
    hold on
    plot(par.time_long(1:end),data.TO,'r^','MarkerSize',8,'linewidth',opt.lw2)
    plot(par.time_long,CO_norm,'bo','MarkerSize',8,'linewidth',opt.lw2)
    area([par.tp,tp_normCT_1],[9999 9999],-9999,'FaceColor','r','FaceAlpha',.1,'EdgeColor','None','linestyle','none');
    xline(par.tp,'k-','linewidth',1.1)
    xline(tp_normCT_1,'k--','linewidth',1.1)
    %legend(par.Tname, [par.Cname,' Normalized'],'Location','Best','interpreter','latex','FontSize',opt.fz1,'box','off')
    xlim([par.tp-1,tp_normCT_1+1]);
    ylim([ymin,ymax])
    set(gca, 'FontSize',opt.fz1);
    ylabel(par.outcome_name ,'interpreter','latex','FontSize',opt.fz2)
    xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
    subplot(2,1,2)
    hold on
    plot(alt_yr,repmat(o_vals.pol_estiCT,[1,sizealt]),'m-','linewidth',opt.lw+1)
    plot(par.time_long,pdiff_CT,'ko','MarkerFaceColor','k','MarkerSize',10)
    plot(tp_normCT_1,gammarCT,'ko','MarkerFaceColor','k','MarkerSize',10)
    plot(par.time_long,zeros(par.nt_long,1),'Color', [0 0 0]+0.2)
    area([par.tp,tp_normCT_1],[9999 9999],-9999,'FaceColor','r','FaceAlpha',.1,'EdgeColor','None','linestyle','none');
    xline(par.tp,'k-','linewidth',1.1)
    xline(tp_normCT_1,'k--','linewidth',1.1)
    %xlim([par.minyear,par.maxyear])
    xlim([par.tp-1,tp_normCT_1+1]);
    ylim([ymin2,ymax2])
    legend('$\gamma$','Location','northwest','FontSize',opt.fz2-1,'interpreter','latex','box','off')
    set(gca, 'FontSize',opt.fz1);
    ylabel('$\gamma(s)$','interpreter','latex','FontSize',opt.fz2)
    xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
    print('-depsc',['output/f8_', par.fapp, '.eps']);
    saveas(f8,['output/f8_', par.fapp, '.png']);
    hold off  

   

end

%% For mapping T to C
%
%% After Normalization T to C 
aux_val = floor(tp_normTC_3);
if tp_normTC_1-par.tp>0
    % its a flip
    aux_val = par.tu;
end
ymin = min([min(data.TO),min(data.CO),min(CO_norm)]);
ymax3 = max([max(TO_norm_TC),max(data.CO)])*1.1;

f9=figure(9);
set(f9, 'Visible',par.visible);
hold on
plot(par.time_long,data.CO,'bo','linewidth',opt.lw)
plot(par.time_long,TO_norm_TC,'r^','linewidth',opt.lw)
plot(par.time_long(1:aux_val),T_norm_TC(1:aux_val),'rx--','linewidth',opt.lw)
plot(par.time_long(1:aux_val),data.C(1:aux_val),'b-','linewidth',opt.lw)
area([tp_normTC_1,par.tp],[9999 9999],-9999,'FaceColor','r','FaceAlpha',.1,'EdgeColor','None','linestyle','none');
xline(par.tp,'k-','linewidth',1.1)
xline(tp_normTC_1,'k--','linewidth',1.1)
legend(par.Cname, [par.Tname,' Normalized'],'Location','Best','interpreter','latex','FontSize',opt.fz2-1,'box','off')
xlim([par.t0,par.t1]);
ylim([ymin,ymax3]);
set(gca, 'FontSize',opt.fz1);
ylabel(par.outcome_name,'interpreter','latex','FontSize',opt.fz2)
xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
print('-depsc',['output/f9_', par.fapp, '.eps']);
saveas(f9,['output/f9_', par.fapp, '.png']);
hold off

% Zoom 
if tp_normTC_1-par.tp>0
    % A flip (Should not be here)
    disp('Flip: You should not be here, check why!')

    f10=figure(10);
    set(f10, 'Visible',par.visible);
    hold on
    plot(par.time_long,data.CO,'bo','MarkerSize',8,'linewidth',opt.lw)
    plot(par.time_long,TO_norm_TC,'r^','MarkerSize',8,'linewidth',opt.lw)
    xline(par.tp,'k-','linewidth',1.1)
    xline(tp_normTC_1,'k--','linewidth',1.1)
    legend(par.Cname, [par.Tname,' Normalized'],'Location','Best','interpreter','latex','FontSize',opt.fz2-1,'box','off')
    set(gca, 'FontSize',opt.fz1);
    ylabel(par.outcome_name ,'interpreter','latex','FontSize',opt.fz2)
    xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
    print('-depsc',['output/f10_', par.fapp, '.eps']);
    saveas(f10,['output/f10_', par.fapp, '.png']);
    hold off    

    f11 = figure(11);
    set(f11, 'Visible',par.visible);
    hold on
    plot(par.time_long,pdiff_TC,'ko','MarkerFaceColor','k','MarkerSize',10)
    plot(par.time_long,zeros(par.nt_long,1),'Color', [0 0 0]+0.2)
    xline(par.tp,'k-','linewidth',1.1)
    xline(tp_normTC_1,'k--','linewidth',1.1)
    legend('$\gamma$','Location','northwest','FontSize',opt.fz2-1,'interpreter','latex','box','off')
    set(gca, 'FontSize',opt.fz1);
    ylabel('$\gamma(s)$','interpreter','latex','FontSize',opt.fz2)
    xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
    print('-depsc',['output/f11_', par.fapp, '.eps']);
    saveas(f11,['output/f11_', par.fapp, '.png']);
    hold off 


else
    
    if floor(tp_normTC_3)<=1
        tp_normTC_3 = floor(tp_normTC_3)+1;
    end
    ymin = min([min(TO_norm_TC(floor(tp_normTC_3)-1:par.tu+1)),min(data.CO(floor(tp_normTC_3)-1:par.tu+1))])*0.95;
    ymax = max([max(TO_norm_TC(floor(tp_normTC_3)-1:par.tu+1)),max(data.CO(floor(tp_normTC_3)-1:par.tu+1))])*1.05;


    f10=figure(10);
    set(f10, 'Visible',par.visible);
    hold on
    plot(par.time_long,data.CO,'bo','MarkerSize',8,'linewidth',opt.lw2)
    plot(par.time_long,TO_norm_TC,'r^','MarkerSize',8,'linewidth',opt.lw2)
    area([tp_normTC_1,par.tp],[9999 9999],-9999,'FaceColor','r','FaceAlpha',.1,'EdgeColor','None','linestyle','none');
    xline(par.tp,'k-','linewidth',1.1)
    xline(tp_normTC_1,'k--','linewidth',1.1)
    legend(par.Cname, [par.Tname,' Normalized'],'Location','Best','interpreter','latex','FontSize',opt.fz2-1,'box','off')
    xlim([tp_normTC_1-1,par.tp+1]);
    ylim([ymin,ymax])
    set(gca, 'FontSize',opt.fz1);
    ylabel(par.outcome_name ,'interpreter','latex','FontSize',opt.fz2)
    xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
    print('-depsc',['output/f10_', par.fapp, '.eps']);
    saveas(f10,['output/f10_', par.fapp, '.png']);
    hold off

    %% Policy Effect Zoom
    alt_yr = [tp_normTC_1,par.tp];
    sizealt=size(alt_yr,2);

    if o_vals.pol_estiCT<0
    ymin2 = min([min(pdiff_TC),o_vals.pol_estiTC])*1.05;
    ymax2 = max([max(pdiff_TC),0.005])*1.05;
        if isnan(ymin2)  
            ymin2 = -0.1; 
        end
    else
    ymin2 = min([min(pdiff_TC),-0.005])*1.05;
    ymax2 = max([max(pdiff_TC),o_vals.pol_estiTC])*1.05;
        if isnan(ymax2)  
            ymax2 = 0.1; 
        end
    end

    f11 = figure(11);
    set(f11, 'Visible',par.visible);
    hold on
    plot(alt_yr,repmat(o_vals.pol_estiTC,[1,sizealt]),'m-','linewidth',opt.lw+1)
    plot(par.time_long,pdiff_TC,'ko','MarkerFaceColor','k','MarkerSize',10)
    plot(tp_normTC_1,(ldr-ldb)/ldb,'ko','MarkerFaceColor','k','MarkerSize',10)
    plot(par.time_long,zeros(par.nt_long,1),'Color', [0 0 0]+0.2)
    area([tp_normTC_1,par.tp],[9999 9999],-9999,'FaceColor','r','FaceAlpha',.1,'EdgeColor','None','linestyle','none');
    xline(par.tp,'k-','linewidth',1.1)
    xline(tp_normTC_1,'k--','linewidth',1.1)
    %xlim([par.minyear,par.maxyear])
    xlim([tp_normTC_1-1,par.tp+1]);
    ylim([ymin2,ymax2])
    legend('$\gamma$','Location','northwest','FontSize',opt.fz2-1,'interpreter','latex','box','off')
    set(gca, 'FontSize',opt.fz1);
    ylabel('$\gamma(s)$','interpreter','latex','FontSize',opt.fz2)
    xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
    print('-depsc',['output/f11_', par.fapp, '.eps']);
    saveas(f11,['output/f11_', par.fapp, '.png']);
    hold off  

    f12 = figure(12);
    set(f12, 'Visible',par.visible);
    subplot(2,1,1)
    hold on
    plot(par.time_long,data.CO,'bo','MarkerSize',8,'linewidth',opt.lw2)
    plot(par.time_long,TO_norm_TC,'r^','MarkerSize',8,'linewidth',opt.lw2)
    area([tp_normTC_1,par.tp],[9999 9999],-9999,'FaceColor','r','FaceAlpha',.1,'EdgeColor','None','linestyle','none');
    xline(par.tp,'k-','linewidth',1.1)
    xline(tp_normTC_1,'k--','linewidth',1.1)
    legend(par.Cname, [par.Tname,' Normalized'],'Location','Best','interpreter','latex','FontSize',opt.fz1,'box','off')
    xlim([tp_normTC_1-1,par.tp+1]);
    ylim([ymin,ymax])
    set(gca, 'FontSize',opt.fz1);
    ylabel(par.outcome_name ,'interpreter','latex','FontSize',opt.fz2)
    xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
    subplot(2,1,2)
    hold on
    plot(alt_yr,repmat(o_vals.pol_estiTC,[1,sizealt]),'m-','linewidth',opt.lw+1)
    plot(par.time_long,pdiff_TC,'ko','MarkerFaceColor','k','MarkerSize',10)
    plot(tp_normTC_1,(ldr-ldb)/ldb,'ko','MarkerFaceColor','k','MarkerSize',10)
    plot(par.time_long,zeros(par.nt_long,1),'Color', [0 0 0]+0.2)
    area([tp_normTC_1,par.tp],[9999 9999],-9999,'FaceColor','r','FaceAlpha',.1,'EdgeColor','None','linestyle','none');
    xline(par.tp,'k-','linewidth',1.1)
    xline(tp_normTC_1,'k--','linewidth',1.1)
    %xlim([par.minyear,par.maxyear])
    xlim([tp_normTC_1-1,par.tp+1]);
    ylim([ymin2,ymax2])
    legend('$\gamma$','Location','northwest','FontSize',opt.fz2-1,'interpreter','latex','box','off')
    set(gca, 'FontSize',opt.fz1);
    ylabel('$\gamma(s)$','interpreter','latex','FontSize',opt.fz2)
    xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
    print('-depsc',['output/f12_', par.fapp, '.eps']);
    saveas(f12,['output/f12_', par.fapp, '.png']);
    hold off   

end % End of flip
%}
end % End of figures
end
%---------------------------------------------------------------------------
function [pphi,x] = make_x(data,opt,m,par)
% Creates Coefficient guesses
    if opt.manual_guess==1
        pphi = m;
        x = pphi;
    elseif opt.manual_guess == 2
         pg = pol_guess(par,data);
        pphi(1,1) = pg(1); % Phi0
        pphi(1,2) = pg(2); % Phi1
        pphi(1,3) = pg(3); % Phi2
        if par.nm ==3
            pphi(1,4) = 0; % Phi3 % Quadratic term
        end
        x = pphi;
    else
      
        % Detect whether series is decreasing or increasing in the basis
        % of the last two points and T
        I = data.T(par.tu)-data.T(1)<0;
        if I ==1 % Debresing
            Cval = data.C(par.tu);
            [~,Tpos] = min(abs(data.T-Cval));% find the closest value in C
            ddl = abs(par.tu-Tpos);
        else % Increasing
            Tval = data.T(par.tu);
            [~,Cpos] = min(abs(data.C-Tval));% find the closest value in C
            ddl = abs(par.tu-Cpos);
        end
        if opt.CT==1
        
        pphi(1,1) = max(data.T)/max(data.C);
        pphi(1,2) = ddl;
        pphi(1,3) = 1;
        else
       
        pphi(1,1) = max(data.C)/max(data.T);
        pphi(1,2) = -ddl;
        pphi(1,3) = 1;
        end
        x = pphi;
    end

    
end
% -------------------------------------------------------------------------
function [fv,pphi,nt] = func_df(x,dta_T,dta_C,opt,par)
% Distance function
[dist,~,~,~,pphi] = func_dist(x,dta_T,dta_C,opt,par);

    nt = length(dist);
    if (opt.wght==1)
        wght = sqrt((1:nt)');
        %wght(nt) = 100.0*wght(nt);
    else
        wght = ones(nt,1);
    end
    wdist = wght.*dist;
  
    fv = 0.5*(wdist')*wdist;    


end
% -------------------------------------------------------------------------
function [dist,pred_C_sc,svec,tvec,pphi] = func_dist(x,dta_T,dta_C,opt,par)
% get the coefficients
pphi = x;

% stage vector corresponding to time before unification
svec = (1:par.tu)';                   % time=stage in T
tvec = func_stage(svec,pphi(:),0,par);    % corresponding time in C

% Predicted scaled values in RO_XX with basic trimming according to C:
% This is the same either from C to T or T to C
   
    I = dta_C>par.cut;
    II = dta_T>par.cut;
    dta_T = dta_T(II);
    aux_svec = svec(II); % save for the graph manubrio
  

try
    pred_C_sc = interp1(tvec(I),dta_C(I),svec(II),'linear',NaN);
    I = isnan(pred_C_sc);
    if sum(I)>1
        %dert_debug =1; 
    end
catch
   
    disp('Interpolation Error')
    pred_C_sc = NaN(size(par.tvec_c(II)));
end

% add the scaling factor
pred_C_sc = pphi(1).*pred_C_sc;


II = isnan(pred_C_sc);
I = not(isnan(pred_C_sc));
if sum(II)>0
    %dert_debug = 1;
end
dist = dta_T(I) - pred_C_sc(I);

if isempty(dist)
    dist = 10e10;
end
% Compute Number of period matched:

if size(dist,1)<par.min_match
    %disp('Not Enough matched points')
    dist = 1e10;
end

% Activar cuando es a manubrio
%{

figure(200)
hold on
plot(aux_svec,pred_C_sc,'b--x')
plot(aux_svec,dta_T,'r--')
plot(svec,dta_C,'b-')
dert_debug = 1;
%}

end
% -------------------------------------------------------------------------
function y = func_stage(x,pphi,back,par)
if back ==0 % in years of C
    y = pphi(2);
    for i = 3:(par.nm+1)
        y = y + (x.^(i-2)).*pphi(i);
    end
else % in years of T
    % Could Compute the inverse, but too comversome
    y = x./pphi(3)- pphi(2);
end

end
%--------------------------------------------------------------------------
function [IDX] = rnd_PS(set,n)

% Draws a random sample from the power set without replacement
% If the number of samples equals to the powerset, then remove the empty
% set
% Learn to remove madrid.
rng(1234)  % set seed
nel = size(set',1);
y = randsample(2^nel,n); 
IDX = de2bi(y);
if size(IDX,2)==nel+1
    % Remove the zeros
    I = IDX(:,end)==1;
    IDX = IDX(:,1:end-1);
    IDX = IDX(not(I),:);

end
vars = 1; 
%function(x) v[intToBits(x) > 0])


end
%--------------------------------------------------------------------------
function cf = spline_reg(x_data,y_data,par,opt)

x_max = max(x_data);
x_min = min(x_data);
if opt.cheb_nodes==1

    [knots] = cheb_roots(par.m);
    x = ((sec(pi/(2*par.m))*knots+1)./2)*(x_max-x_min)+x_min;
    knots = x;
elseif opt.cheb_nodes ==2
    % Evenly spaced 
    x = linspace(x_min,x_max,par.m);
    knots = augknt(x,4);
else
    knots = par.m; % Automatic by matlab
end
% 4 to have 2 continous derivatives
cf_func = spap2(knots,par.n,x_data,y_data);
cf = fnval(cf_func,x_data);


end
%--------------------------------------------------------------------------
function [theta_guess,fx] = guess_theta(x_data,y_data,par,opt)
% To find a guess for theta I do the following: 
% I fit a polinomial to the data, with this I get continous values. 
% Then I use thos values to solve for analytical thetas which are my guess:
% 

x_max = max(x_data);
x_min = min(x_data);
zi = cheb_roots(par.m);
% Step 5: Convert nods to x interval x_max,x_min
knots = (sec(pi/(2*par.m))*zi+1)/2*(x_max-x_min)+x_min;
    
fx = csaps(x_data,y_data,par.sp,knots);
%fx_func = spap2(par.m,par.n+10,x_data,y_data);
%fx = fnval(fx_func,knots);
theta_guess = NaN(par.n,1);
val_aux = NaN(par.m,1);
% Here I compute the analitical thetas accordinf to Nakajima:
Taux = Tox(par.n-1,zi);
for j = 1:par.n
    if j ==1
        theta_guess(j) = 1/par.m*(sum(fx));
    else
        for i = 1:par.m
            %val_aux(i) = fx(i).*chebyshevT(j-1,zi(i));
            %val_aux(i) = fx(i).*Tox(j-1,zi(i));
            val_aux(i) = fx(i).*Taux(i,j);
        end
        I = isnan(val_aux);
        if sum(I)>0
            disp('some error on polinomial')
        end
        theta_guess(j) = (2/par.m)*sum(val_aux);
    end
end



%{
dert_debug = 1;
figure(101)
hold on
plot(x_data,y_data,'k-o','markerfacecolor','k','linewidth',0.5)
plot(knots,fx,'b-','linewidth',1.1)
%plot(knots,gx,'rx-','linewidth',1.1)
%}
end
%--------------------------------------------------------------
function T = Tox(or,x)
%{
Manual Chabyshev polinomials, cuz its faster:
%}
%{
    if or ==0
        val = 1;
    elseif or ==1
        val = x;
    elseif or ==2
        val = 2*(x^2)-1;
    elseif or ==3
        val = 4*(x^3)-3*(x);
    elseif or ==4
        val = 8*(x^4)-8*(x^2)+1;
    elseif or ==5
        val = 16*(x^5)-20*(x^3)+5*(x);
    elseif or ==6
        val = 32*(x^6)-48*(x^4)+18*(x^2)-1;
    elseif or ==7      
        val = 64*(x^7)-112*(x^5)+56*(x^3)-7*x;
    elseif or ==8      
        val = 128*(x^8)-256*(x^6)+160*(x^4)-32*(x^2)+1;
    end
    %}
n = or+1;    
m = length(x);
T = ones(m, 1);
dT = zeros(m,1);
ddT = zeros(m,1);
for i = 2 : n
    if i == 2
        T = [T, x];
        dT = [dT, ones(m,1)];
        ddT = [ddT, zeros(m,1)];
    else
        T = [T, 2.0 * x .* T(:, i-1) - T(:, i-2)];
        dT = [dT, 2.0 * T(:, i-1) + 2.0 .* x.* dT(:, i-1) - dT(:, i-2)];
        ddT = [ddT, 2.0 * dT(:, i-1) + 2.0 .* dT(:, i-1) + 2.0 .* x.* ddT(:, i-1) - ddT(:, i-2)];
    end
end

end
%--------------------------------------------------------------------------
function [roots_vec] = cheb_roots(m)
%{
Computes the roots of the chebychev polinomial of order m
%}
i_aux = (1:m)';
roots_vec = -cos((2.*i_aux-1).*pi./(2*m));

end
%--------------------------------------------------------------------------
function [fv,chebyval] =  dist_cheby(theta,x_data,y_data,par)
% Evaluate the Cheby polinomial


% Step 1: 
x_max = max(x_data);
x_min = min(x_data);
% Step 2: Pick the order of the Chebyshev Polinomial
n = par.n;
nobs = size(x_data,1);
% Step 3: Pick the number of Collocation points
m = par.m;  %par.nobs;% Number of collocation points same as order or more

% Compute solution thetas, I need data counterpart to do this: 
% So I just evaluate the Polinomial 
% Convert data points x to -1,1 
zi = (1/sec(pi/(2*m))).*(((2*(x_data-x_min))/(x_max-x_min))-1);

cheby_poli = NaN(n,nobs);
Taux = Tox(n-1,zi);
for i = 1:n
    for ii = 1:nobs
        %cheby_poli(i,ii) = theta(i)*chebyshevT(i-1,zi(ii));
        %cheby_poli(i,ii) = theta(i)*Tox(i-1,zi(ii));
        cheby_poli(i,ii) = theta(i)*Taux(ii,i);
    end
end

%z = ((sec(pi/(2*par.m))*roots_vec+1)./2)*(z_max-z_min)+z_min;
%z = (1/sec(pi/(2*m))).*(((2*(x_data-x_min))/(x_max-x_min))-1);

I = isnan(cheby_poli);
if sum(I,'all')>0
    disp('Something wrong is happening')
end
chebyval = sum(cheby_poli,1)';
% Get the difference: 
I = not(isnan(chebyval));
dist = y_data(I) - chebyval(I);

if isempty(dist)
    disp('Why is dist empty')
    dist = 10e10;
end

 wdist = dist;
  
 fv = 0.5*(wdist')*wdist;   

end
%--------------------------------------------------------------------
function [cf,theta_sol,fval]=coll(x_data,y_data,par,opt)
%{
x, y = data input, does not allow for missing values
n = degree of the polinomial n_pol-1 = order of the polinomial

%}


% Guess theta (Cheby Coeffs)

[theta_x0,~] = guess_theta(x_data,y_data,par,opt);



[theta_sol,fval,ef] =  fminunc(@(theta) dist_cheby(theta,x_data,y_data,par),theta_x0,opt.optmin);

[~,cf] = dist_cheby(theta_sol,x_data,y_data,par);

end
%--------------------------------------------------------------------------
function pg = pol_guess(par,data)
% This fits a polinomial prepolicy and uses the analitical coefficients as
% guesses
% Chooses the guess with positive phi(2)

    pDC = polyfit(par.time(1:par.tu),data.C(1:par.tu),2);
    pDT = polyfit(par.time(1:par.tu),data.T(1:par.tu),2);
    auxDC = polyval(pDC,par.time(1:par.tu));
    auxDT = polyval(pDT,par.time(1:par.tu));
    % Normalize analically
    pa.A = (pDT(1)*pDC(2))/(pDC(1)*pDT(2));
    pa.B = (pDT(1)*2)/pDT(2);
    pa.C = pDT(2)/pDT(3);

    pa.c1 = pDC(2)*pa.A-pDC(3)*pa.C;
    pa.c2 = pDC(2)*pa.B+2*pDC(1)*pa.A-pDC(2)*pa.C;
    pa.c3 = 2*pDC(1)*pa.B-pDC(1)*pa.C;

    rvals = roots([pa.c3,pa.c2,pa.c1]);

    pa.phi1 = rvals(2); % choose this because it gives a positive phi2
    pa.phi2 = pa.A+pa.B*pa.phi1;
    pa.phi0 = (pa.phi2^2)*pDC(1)/pDT(1);

    if pa.phi2<0
        pa.phi1 = rvals(1); % choose this because it gives a positive phi2
        pa.phi2 = pa.A+pa.B*pa.phi1;
        pa.phi0 = (pa.phi2^2)*pDC(1)/pDT(1);
    end

    if pa.phi2<0
        %if opt.CT==1
        pa.phi0 = 1.1;
        pa.phi1 = -6; % choose this because it gives a positive phi2
        pa.phi2 = 1;

        %else
        %pa.phi1 = 6; % choose this because it gives a positive phi2
        %pa.phi2 = 1;
        %pa.phi0 = 0.9;
        %end
        if par.iter==0
        disp('Polinomial guess failed: using default guess')
        end
    end
    pg(1) = 1/pa.phi0;
    pg(2) = -pa.phi1/pa.phi2;
    pg(3) = 1/pa.phi2;
    
    if abs(pg(2))>2*par.nt
        %if opt.CT==1
        pg(1) = 0.8;
        pg(2) = 2; % choose this because it gives a positive phi2
        pg(3) = 1;

        %else
        %pg(2) = -6; % choose this because it gives a positive phi2
        %pg(3) = 1;
        %pg(1) = 0.8;
        %end
        if par.iter==0
        disp('Using default guess instead of polinomial guess')
        end
    end

   
    %{
    dert_debug = 1;
    tnorm = pg(2)+pg(3)*par.time(1:par.tu);
    figure(421)
    hold on
    plot(par.time(1:par.tu),auxDC,'b-')
    plot(par.time(1:par.tu),auxDT,'ro-')
    plot(tnorm,auxDC.*pg(1),'bx-')
    %

    %
    figure(401)
    hold on
    plot(par.time,data.C,'b--')
    plot(par.time,data.T,'r--')
    plot(par.time(1:par.tu),auxDC,'bx-')
    plot(par.time(1:par.tu),auxDT,'rx-')
    %ylim([-50,250])
    %xlim([0,80])
    %}

end
%------------------------------------------------------------------------
function [S_mean,S_median,UCI,BCI,ntest,nobs,S_Edist]=getCI(ser,lev,flag)
%{
This function computes the confidence bands for a given empirical
distribution using the percentile method.

This function also performs the loops
nc: number of Columns, 
nb: number of Boorstrap iterations
ser: series, should be NBxNC
flag: if 1 it is an outlier and should be taken out
nobs: number of healthy iterations-NaN's
%}
par.nc = size(ser,2);
%par.nb = nb;
S_Edist = [];              % Empirical distribution for histogram
S_median = NaN(par.nc,1);
S_mean   = NaN(par.nc,1);
UCI = NaN(par.nc,1);
BCI = NaN(par.nc,1);
CI_BU = NaN(par.nc,2);     % (B)Bottom (U)Upper
ntest = NaN(par.nc,1);     % Jb test
nobs = NaN(par.nc,1);

I = logical(flag);
ser = ser(not(I),:);

for rc = 1:par.nc
    sort_gam = sort(ser(:,rc));
    % exclude the NaNs, if only NaN's then skip
    I = isnan(sort_gam);
    if sum(I)==size(sort_gam,1)
        continue
    end
    sort_gam = sort_gam(not(I));
    S_Edist = sort_gam;
    nb_t = size(sort_gam,1);
    nobs(rc,1) = nb_t;   % Number of observations, this excludes the nans
    
    cut = round(nobs(rc,1)*lev/2); % cut for computing confidence bounds
    if cut ==0
        %{
        CI_BU(rc,:) = [-1,1];
        UCI(rc,1) = 1;
        BCI(rc,1) = -1;
        S_median(rc,1) = 0;
        S_mean(rc,1) = 0;
        %}
        CI_BU(rc,:) = [NaN,NaN];
        UCI(rc,1) = NaN;
        BCI(rc,1) = NaN;
        S_median(rc,1) = NaN;
        S_mean(rc,1) = NaN;
    else
        CI_BU(rc,:) = [sort_gam(cut),sort_gam(nb_t-cut+1)];
        UCI(rc,1) = CI_BU(rc,2);
        BCI(rc,1) = CI_BU(rc,1);
        S_median(rc,1) = median(sort_gam);
        S_mean(rc,1) = mean(sort_gam);
        try
           % JB test with H0: Normal
            %[~,p_jb] = jbtest(sort_gam);
            %ntest(rc,1) = p_jb;%nb_t; 
        catch
        end
    end

end


end
%--------------------------------------------------------------------------
function [Edist_peA,Edist_pe,Edist_tpn,med_tpn,med_peA,y_min,y_max] = bb_perform(datas,data,out_r,par,opt)

    disp('Initializing Boostrap')
    opt.warn = 0;
    if opt.warn ==0
        disp('Warnings and iter counter turned off for boostrap iters')
    end
    %% Back out Residuals for Bootstrapping
    
    out_r.Csmooth = datas.C;
    out_r.Tsmooth = datas.T;
    out_r.resid_C = zeros(par.nt,1); % for when there are zero values
    out_r.resid_T = zeros(par.nt,1); % for when there are zero values
    I =  datas.CO>0;
    J =  datas.C>0;
    II = logical(I.*J);
    %II = 1:par.nt>10; % Eliminate Residual Outliers
    out_r.resid_C(II) = log(datas.CO(II))-log(datas.C(II));
    I =  datas.TO>0;
    J =  datas.T>0;
    II = logical(I.*J);
    %II = 1:par.nt>10;
    out_r.resid_T(II) = log(datas.TO(II))-log(datas.T(II));
    if opt.smooth==1
        out_r.resid_C = out_r.resid_C(1:par.tu);
        out_r.resid_T = out_r.resid_T(1:par.tu);
    end
 
    
    % Fish for >0 autocorrelation in C
    varC = var(out_r.resid_C);
    rhoC = corr(out_r.resid_C(2:end),out_r.resid_C(1:end-1));
    % Fish for >0 autocorrelation in T
    varT = var(out_r.resid_T);
    rhoT = corr(out_r.resid_T(2:end),out_r.resid_T(1:end-1));
    if rhoT>0 || rhoC>0
        bb_d = 1;
        disp('Using Block Bootrap...')
        disp(['Residual Autocorr C: ',num2str(rhoC)])
        disp(['Residual Autocorr T: ',num2str(rhoT)])
    else
        bb_d = 0;
        disp(['Residual Autocorr C: ',num2str(rhoC)])
        disp(['Residual Autocorr T: ',num2str(rhoT)])
        %disp('Using Block Bootrap...')
    end
 
    
    I = out_r.resid_C ==0;
    II = out_r.resid_T ==0;    
    [~,t0A] = max(not(I));
    [~,t0B] = max(not(II));
    ti = max([t0A,t0B]); % Select the region which zeros start later
    
   
    
    %}
  %% Initialize for Bootstrapping:
    flag_out = zeros(par.nb,1);                % Outlier indicator 1: outlier, 0: healthy iteration   
    flag_out_sign = zeros(par.nb,1);           % Outlier indicator 1: outlier, 0: healthy iteration
    flag_out_flip = zeros(par.nb,1);           % Outlier indicator 1: outlier, 0: healthy iteration
    flag_out_nons = zeros(par.nb,1);           % Outlier indicator 1: outlier, 0: healthy iteration
    phi_bd = -99999*ones(par.nb,par.nm+1);            % Coefficients
    olp_bd = -99999*ones(par.nb,1);            % Identification interval
    pol_estiCT = -99999*ones(par.nb,1);        % Gamma
    LS_CT = -99999*ones(par.nb,1);             % LS
    fvalCT = -99999*ones(par.nb,1);            % Norm Min Error
    pdiff_CT = -99999*ones(par.nt,par.nb);      
    pdiff_CT_int = -99999*ones(par.nt,par.nb); % T Interpolated
    pdiff_CT_pre = -99999*ones(par.nt,par.nb); % 
    pdiff_CT_int_pre = -99999*ones(par.nt,par.nb); % Prepolicy T Interpolated
    CO = -99999*ones(par.nt,par.nb);           % C Series
    CO_norm = -99999*ones(par.nt,par.nb);      % C Series Normalized
    TO = -99999*ones(par.nt,par.nb);           % T Series
    time_d = -99999*ones(par.nt,par.nb);       % Time
    time_d_norm = -99999*ones(par.nt,par.nb);  % Time normalized
    tp_normCT = -99999*ones(par.nb,1);         % Normalized policy date
    
    T = NaN(par.nt,par.nb);                 % T Smooth Series
    C = NaN(par.nt,par.nb);                 % C Smooth Series
    Cnorm = NaN(par.nt,par.nb);             % C Smooth Series
    s = RandStream('mlfg6331_64');
    rng(123,'twister')
    %
    if opt.smooth==0
        opt.b=0;
    end
    if opt.b==1
        disp('Boostrapping...')
        par.iter = 1; % discard first iteration for display purposes
     
        for i = 1:par.nb
          
            par.vis_ind = 0; % Turn off graphs
            close all
            cmeC_bd = zeros(par.nt,1);    % Sampled residuals C
            cmeT_bd = zeros(par.nt,1);    % Sampled residuals T
            cmeC_bd_bb = zeros(par.nt,1); % Sampled residuals for Block Boostrap
            cmeT_bd_rb = zeros(par.nt,1); % Sampled residuals for Traditional Boostrap


            rv_draw_pre = ones(par.nt,1); % Draw Indicator Pre_polici
            rv_draw = ones(par.nt,1);     % Draw Indicator Post_policy only if smooth ==2

            %t0 = 1; % The region that stards further. 
            nr = par.tu-ti+1;             % Number of prepolicy days
            nx = ceil(nr/par.wdth);       % Number of bootstrap intervals in prepolicy days ex 5
            ind = randi([0 nr-par.wdth],nx,1); % Generate nx=5 draws
            tvec = [ti:par.tu]';
            for xc=1:nx                   % Roll on boothstrap intervals
                for tc=1:par.wdth %       Roll on values inside an interval
                    tcc = (xc-1)*par.wdth+tc+(ti-1); % 
                    if ( tcc>par.tu )
                        break
                    end
                    rv_draw_pre(tcc) = tvec(ind(xc)+tc);
                end
            end
            % after policy reform:
            if opt.smooth==2
                % Draw only from prepolicy error 
                nr2 = par.nt-(par.tu+1)+1; % number of post policy days 
                nx = ceil(nr2/par.wdth);
                ind = randi([0 nr-par.wdth],nx,1);
                tvec = [ti:par.tu]';
                for xc=1:nx
                    for tc=1:par.wdth
                        tcc = (xc-1)*par.wdth+tc+par.tu;
                        if ( tcc>par.nt )
                            break
                        end
                        rv_draw(tcc) = tvec(ind(xc)+tc);
                    end
                end

            end


            % Pre-policy
            % Form Block Bootstrap
            cmeC_bd_bb(ti:par.tu) =  out_r.resid_C(rv_draw_pre(ti:par.tu));
            cmeT_bd_bb(ti:par.tu) =  out_r.resid_T(rv_draw_pre(ti:par.tu));

            % Sample with replacement
            cmeC_bd_rb(ti:par.tu) = (randsample(s,out_r.resid_C(ti:par.tu)',par.tu-ti+1,true))';
            cmeT_bd_rb(ti:par.tu) = (randsample(s,out_r.resid_C(ti:par.tu)',par.tu-ti+1,true))';

            if bb_d==1 % Block bootstrap 
                cmeC_bd(ti:par.tu) = cmeC_bd_bb(ti:par.tu);
                cmeT_bd(ti:par.tu) = cmeT_bd_bb(ti:par.tu);
            else
                cmeC_bd(ti:par.tu) = cmeC_bd_rb(ti:par.tu);
                cmeT_bd(ti:par.tu) = cmeT_bd_rb(ti:par.tu);
            end




             % Post-Policy

            cmeC_bd(par.tu+1:par.nt) = 0;%out_r.resid_C(rv_draw(par.tp+1:par.nt));
            cmeT_bd(par.tu+1:par.nt) = 0;%out_r.resid_T(rv_draw(par.tp+1:par.nt)); 

            data.C = out_r.Csmooth.*exp(cmeC_bd);%  Dashed Blue    
            data.T = out_r.Tsmooth.*exp(cmeT_bd);%  Dashed Red

            % Smooth 
            [datas]=smooth_ts(data,par,opt);
            %% Perform SBI

                data.T = datas.T;
                data.C = datas.C;
                data.TO = datas.TO;
                data.CO = datas.CO;


             % Graph Smooth
             %{
                if opt.graph_smooth==1
                    figure(1000)
                    hold on
                    plot(par.time,datas.T,'r-')
                    plot(par.time,datas.C,'b-')
                    plot(par.time,datas.TO,'rx')
                    plot(par.time,datas.CO,'bo')
                end
            %}

            [out] = norm_func(data,par,opt);
        %% Trim the sample: 
        %{
         ut.flag_out_nons ==1 means the algorithm is finding nonsence
         Nonsence can be:(1)extreme cases (2)non-convergence (3)gamma = NaN
         I want to trim out nonsence and outliers
         Outliers are:   
         -crazy psis (detected by default as an outlier)
        %}
                % trim nonsence
                if isnan(out.pol_estiCT)
                    flag_out(i) = 1;
                end
                if out.flag_out_nons ==1  
                    flag_out(i) = 1;
                end
                % Dont keep the flips
                if out.flag_out_flip ==1
                    flag_out(i) = 1;
                end
                if out.flag_out_sign ==1
                    %flag_out(i) = 1;
                end
                flag_out_nons(i) = out.flag_out_nons;
                flag_out_sign(i) = out.flag_out_sign;
                flag_out_flip(i) = out.flag_out_flip;

            if flag_out(i)==0
                par.iter = par.iter+1;
                if par.iter>(par.nb_lim+1) % 500 bootstrap iters
                    flag_out(i:end) = 1;
                    break
                end
                if opt.warn==1
                disp(['Bootstrap iteration No: ',num2str(par.iter-1)])  
                end
                phi_bd(i,:) = out.psi_estiCT;           % Coefficients
                olp_bd(i) = out.olpCT;                  % Overlap interval
                pol_estiCT(i) = out.pol_estiCT;         % Gamma
                LS_CT(i) = out.LS_CT;                   % LS
                fvalCT(i) = out.fvalCT;                 % Norm Min Error
                pdiff_CT(:,i) = out.pdiff_CT;           % Cum Gamma
                pdiff_CT_int(:,i) = out.pdiff_CT_int;   % Cum Gamma Interpolated
                pdiff_CT_pre(:,i) = out.pdiff_CT_pre;   % Cum Gamma PreP 
                pdiff_CT_int_pre(:,i) = out.pdiff_CT_int_pre; % Cum Gamma PreP Interpolated
                CO(:,i) = out.CO;                 % C Series
                CO_norm(:,i) = out.CO_norm;       % C Series Normalized
                TO(:,i) = out.TO;                 % T Series
                time_d(:,i) = out.time;           % Time
                time_d_norm(:,i) = out.time_norm; % Time Normalized
                tp_normCT(i) = out.tp_normCT;     % Norm policy date

                auxsN = size(out.C,1); 
                C(1:auxsN,i) = out.C;               % C Smooth Series
                auxsN = size(out.Cnorm,1);
                Cnorm(1:auxsN,i) = out.Cnorm;       % C Smooth Series Normalized
                auxsN = size(out.T,1);
                T(1:auxsN,i) = out.T;               % T Smooth Series 
            end

        end % end boostrap iter
        if par.iter<(par.nb_lim+1)
            warning(['Only ',num2str(par.iter-1),' out of ',num2str(par.nb_lim),' converged'])
        end
    else
        if opt.smooth==0
            disp('No smoothing step: Unable to boostrap')
        end
        if opt.smooth==1 && opt.b==0
            disp('User chose to skip boostrap step')
        end
        flag_out = ones(par.nb,1);                % Outlier indicator 1: outlier, 0: healthy iteration  
    end % end of boostrap conditional
    %% Get CI's throught the percentile method
   
    %% See how many outliers I am eliminating
    Io1 = flag_out_flip ==1;
    Io2 = flag_out_sign ==1;
    Io3 = flag_out_nons ==0;
    I1A = Io1.*Io2.*Io3 ;
    I1B = Io2.*Io3;
    no_ex = nansum(I1A);
    rat_ex = nansum(I1A)/nansum(I1B);
    %% Get the CI's
    %par.lev = 0.1; % 90% 
    % Phi's
    [mean_phi,med_phi,UCI_phi,BCI_phi,~,nobs_phi,Edist_phi]=getCI(phi_bd,par.lev,flag_out);
    % Policy Effects
    [mean_pe,med_pe,UCI_pe,BCI_pe,~,nobs_pe,Edist_pe]=getCI(pol_estiCT,par.lev,flag_out);
       % Cum gamma
    [mean_cg,med_cg,UCI_cg,BCI_cg,~,nobs_cg,Edist_cg]=getCI(pdiff_CT_int',par.lev,flag_out);
     % Series 
    [mean_cd,med_cd,UCI_cd,BCI_cd,~,nobs_cd,Edist_cd]=getCI(CO',par.lev,flag_out);
    [mean_td,med_td,UCI_td,BCI_td,~,nobs_td,Edist_td]=getCI(TO',par.lev,flag_out);
    % policy date
    [mean_tpn,med_tpn,UCI_tpn,BCI_tpn,ntest_tpn,nobs_tpn,Edist_tpn]=getCI(tp_normCT,par.lev,flag_out);
    
    flag_out_altA = flag_out;
    flag_out_altB = flag_out;
    % Cut around the median Estimate idenfication window
    BI = (tp_normCT-par.tp)>(med_tpn-par.tp)+(med_tpn-par.tp).*par.wd;
    BII = (tp_normCT-par.tp)<(med_tpn-par.tp)-(med_tpn-par.tp).*par.wd;
    % Cut for day around the Point Estimate Identification window  
    AI = (tp_normCT-par.tp)>((out_r.tp_normCT-par.tp)+(out_r.tp_normCT-par.tp).*par.wd);
    AII = (tp_normCT-par.tp)<((out_r.tp_normCT-par.tp)-(out_r.tp_normCT-par.tp).*par.wd);

    
    flag_out_altA(AI) = 1;
    flag_out_altA(AII) = 1;
  
    [mean_peA,med_peA,UCI_peA,BCI_peA,~,nobs_peA,Edist_peA]=getCI(pol_estiCT,par.lev,flag_out_altA);
    
    % Get it for the diff graph:
    [mean_cgA,med_cgA,UCI_cgA,BCI_cgA,~,~,Edist_cgA]      = getCI(pdiff_CT',par.lev,flag_out_altA);     % 
    [mean_cgPA,med_cgPA,UCI_cgPA,BCI_cgPA,~,~,Edist_cgPA] = getCI(pdiff_CT_pre',par.lev,flag_out_altA); % 
 
    flag_out_altB(BI) = 1;
    flag_out_altB(BII) = 1;
    
    [mean_peB,med_peB,UCI_peB,BCI_peB,~,nobs_peB,Edist_peB]=getCI(pol_estiCT,par.lev,flag_out_altB);
    
    % Get it for the diff graph:
    [mean_cgB,med_cgB,UCI_cgB,BCI_cgB,~,~,Edist_cgB]      = getCI(pdiff_CT',par.lev,flag_out_altB);     % 
    [mean_cgPB,med_cgPB,UCI_cgPB,BCI_cgPB,~,~,Edist_cgPB] = getCI(pdiff_CT_pre',par.lev,flag_out_altB); % 
 
  
    % Get the case bordering the median effect
    IC = (tp_normCT-par.tp) > (med_tpn-par.tp)-0.05 & (tp_normCT-par.tp) <(med_tpn-par.tp)+0.05;
    [~,pos_med] = max(IC);
    %}
    %% Save results 
       save(['output/',par.name_fapp,par.fapp]) 
       %load(['output/',par.name_fapp,par.fapp]) 
    %% Results table and save it
    stat = {'gamma bench';'gamma ard bench window';'gamma ard median window ';'windon size'};
    for ii = 0:par.nm
        stat = [stat;{['phi_',num2str(ii)]}];
    end
   estimate = [out_r.pol_estiCT;out_r.pol_estiCT;mean_peB;out_r.olpCT;out_r.psi_estiCT(:)];
   mean   = [mean_pe;mean_peA;mean_peB;mean_tpn-par.tp;mean_phi(:)];
   median = [med_pe;med_peA;med_peB;med_tpn-par.tp;med_phi(:)];
   lower_bound  = [BCI_pe;BCI_peA;BCI_peB;BCI_tpn-par.tp;BCI_phi(:)];
   upper_bound  = [UCI_pe;UCI_peA;UCI_peB;UCI_tpn-par.tp;UCI_phi(:)];
   nobs   = [nobs_pe;nobs_peA;nobs_peB;nobs_tpn;nobs_phi(:)];
   
   table1 = table(stat,estimate,mean,median,lower_bound,upper_bound,nobs);
   disp(['Computing ',num2str((1-par.lev)*100),' CIs'])
   disp(table1)
    save('output/results_table','table1')
    if opt.b==1
    %% Figures
        if max(med_cgA)>0
            y_max = max(max([max(UCI_cgA),max(UCI_cgPA(par.tu-3:par.tu))])).*1.05;
            y_min = min([min(BCI_cgA),min(BCI_cgPA(par.tu-3:par.tu)),0]).*1.05;
        else
            y_max = max([max(UCI_cgA),max(UCI_cgPA(par.tu-3:par.tu)),0]).*1.05;
            y_min = min(min([min(BCI_cgA),min(BCI_cgPA(par.tu-3:par.tu))])).*1.05;
        end
         mean_peB = interp1([floor(out_r.tp_normCT-par.t0+1),ceil(out_r.tp_normCT-par.t0+1)],[med_cgA(floor(out_r.tp_normCT-par.t0+1)),med_cgA(ceil(out_r.tp_normCT-par.t0+1))],out_r.tp_normCT-par.t0+1,'linear',NaN);
        if isnan(mean_peB)
             mean_peB = med_cgA(floor(out_r.tp_normCT-par.t0+1));
        end
    

    %% Cumulative gamma
   
    figure(20);
    hold on
    grid on
    plot([par.tp,out_r.tp_normCT],[mean_peB,mean_peB],'m-','linewidth',opt.lw2)
    plot(out_r.tp_normCT,mean_peB,'o','Color','m','MarkerSize',6,'linewidth',2.2);
    area([par.tp,out_r.tp_normCT],[9999 9999],-9999,'FaceColor','r','FaceAlpha',.1,'EdgeColor','None','linestyle','none');
    p1 = plot(par.t0:par.t1,med_cgA,'m-','linewidth',opt.lw2);
    p2 = plot(par.t0:par.t1,mean_cgA,'m--','linewidth',opt.lw2);
    p4 = plot(par.t0:par.t1,UCI_cgA,'k--','linewidth',opt.lw);
    plot(par.t0:par.t1,BCI_cgA,'k--','linewidth',opt.lw)
    %p5 = plot(out_data.time_norm,out_data.pdiff_CT_int,'o','MarkerSize',8,'linewidth',2.2,'Color',[0.4940 0.1840 0.5560]);
    %p6 = plot(par.time,out_data.pdiff_CT,'^','MarkerSize',8,'linewidth',2.2,'Color',[0.4940 0.1840 0.5560]);
    %plot(out_data.time_norm,out_data.pdiff_CT_int_pre,'ko','MarkerSize',8,'linewidth',2.2);
    %plot(par.time,out_data.pdiff_CT_pre,'k^','MarkerSize',8,'linewidth',2.2);
    plot(par.t0:par.t1,mean_cgPA,'k-','linewidth',opt.lw2);
    plot(par.t0:par.t1,UCI_cgPA,'k--','linewidth',opt.lw);
    plot(par.t0:par.t1,BCI_cgPA,'k--','linewidth',opt.lw);
    % Connect all dots
    plot([par.tp,par.tp+1],[mean_cgPA(par.tu),mean_cgA(par.tu+1)],'m--','linewidth',opt.lw2) % Center
    plot([par.tp,par.tp+1],[mean_cgPA(par.tu),med_cgA(par.tu+1)],'m-','linewidth',opt.lw2) % Center
    plot([par.tp,par.tp+1],[UCI_cgPA(par.tu),UCI_cgA(par.tu+1)],'k--','linewidth',opt.lw) % Center
    plot([par.tp,par.tp+1],[BCI_cgPA(par.tu),BCI_cgA(par.tu+1)],'k--','linewidth',opt.lw) % Center
    ylim([-0.11,y_max])
    %xlim([par.tp-3,med_tpn])
    xlim([par.tp-3,out_r.tp_normCT+1])
    yline(0,'-','Color',[0.5,0.5,0.5],'linewidth',opt.lw)
    xline(par.tp,'-','Color',[0.5,0.5,0.5],'linewidth',opt.lw)
    xline(out_r.tp_normCT,'--','Color',[0.5,0.5,0.5],'linewidth',opt.lw)
    %xline(med_tpn,'--','Color',[0.5,0.5,0.5],'linewidth',opt.lw)
    %legend([p5,p1,p2,p4],{'Point Estimate','Median','Mean','90% CI'},'location','best')
    legend([p1,p2,p4],{'Median','Mean','$90\%$ CI'},'location','best','FontSize',opt.fz2,'interpreter','latex','box','off')
    set(gca, 'FontSize',opt.fz1);
    %title('Policy Effect ($\gamma_{s}$)','interpreter','latex')
    ylabel('$\gamma(s)$','interpreter','latex','FontSize',opt.fz2)
    xlabel('Stage','interpreter','latex','FontSize',opt.fz2)

    %% The Regions 
    
    data.rate(:,1) = data.oriC; %  Dashed Blue    
    data.rate(:,2) = data.oriT;
    figure(21);
    hold on
    set(gca, 'FontSize',opt.fz1);
    p1 = plot(par.t0:par.t1,data.rate(:,1),'bo','linewidth',1.1);
    p2 = plot(par.t0:par.t1,data.rate(:,2),'r^','linewidth',1.1);
    p3 = plot(par.t0:par.t1,med_cd,'b-','linewidth',1.1);
    p4 = plot(par.t0:par.t1,med_td,'r-','linewidth',1.1);
    p5 = plot(par.t0:par.t1,UCI_cd,'b--','linewidth',1.1);
    p6 = plot(par.t0:par.t1,BCI_cd,'b--','linewidth',1.1);
    p7 = plot(par.t0:par.t1,UCI_td,'r--','linewidth',1.1);
    p8 = plot(par.t0:par.t1,BCI_td,'r--','linewidth',1.1);
    p9 = plot([1,2],[NaN,NaN],'k--','linewidth',1.1);
    xline(par.tp,'k-','linewidth',1.1);
    %ylim([0,ym])
    xlim([par.t0,par.t1])
    legend([p3,p4,p9],{'Median $\mathcal{C}$','Median $\mathcal{T}$','$90\%$ CI'},'location','best','FontSize',opt.fz2,'interpreter','latex','box','off')
    xlabel(par.time_name,'interpreter','latex','FontSize',opt.fz2)
    ylabel(par.outcome_name,'interpreter','latex','FontSize',opt.fz2)      
        
    
    figure(22)
    hold on
    set(gca, 'FontSize',opt.fz1);
    p1 = plot(par.time(1:par.tu),CO(1:par.tu,pos_med),'bo','linewidth',1);
    p2 = plot(par.time(1:par.tu),TO(1:par.tu,pos_med),'r^','linewidth',1);
    %plot(par.tp+1:par.nt,CO(par.tp+1:par.nt,pos_med),'o','linewidth',0.7,'Color',[0,0,1,0.4]);
    %plot(par.tp+1:par.nt,TO(par.tp+1:par.nt,pos_med),'^','linewidth',0.7,'Color',[1,0,0,0.4]);
    p12 = plot(par.time,C(:,pos_med),'b-','linewidth',1.4);
    p13 = plot(par.time,T(:,pos_med),'r-','linewidth',1.4);
    %p10 = plot(1:par.nt,dat.dC(:,1,2),'b-','linewidth',1.1,'Color',[0,0,1,0.4]);
    %p11 = plot(1:par.nt,dat.dT(:,1,2),'r-','linewidth',1.1,'Color',[1,0,0,0.4]);
    %p3 = plot(1:par.nt,med_cd,'b-','linewidth',1.1);
    %p4 = plot(1:par.nt,med_td,'r-','linewidth',1.1);
    p5 = plot(par.t0:par.t1,UCI_cd,'--','linewidth',1,'Color',[0,0,1,0.5]);
    p6 = plot(par.t0:par.t1,BCI_cd,'--','linewidth',1,'Color',[0,0,1,0.5]);
    p7 = plot(par.t0:par.t1,UCI_td,'--','linewidth',1,'Color',[1,0,0,0.5]);
    p8 = plot(par.t0:par.t1,BCI_td,'--','linewidth',1,'Color',[1,0,0,0.5]);
    p9 = plot([1,2],[NaN,NaN],'k--','linewidth',1.1);
    xline(par.tp,'k-','linewidth',1.1);
    %ylim([0,ym])
    xlim([par.t0,par.t1])
    legend([p12,p13,p9],{'Median $\mathcal{C}$','Median $\mathcal{T}$','$90\%$ CI'},'location','best','FontSize',opt.fz2,'interpreter','latex','box','off')
    xlabel(par.time_name,'interpreter','latex','FontSize',opt.fz2)
    ylabel(par.outcome_name,'interpreter','latex','FontSize',opt.fz2)
    
    %% Restricted 
    
    % Around the benchmark window 
    figure(23);
    hold on
    set(gca, 'FontSize',opt.fz1);
    histogram(Edist_peA,'BinWidth',0.1/2.5,'FaceColor',[0.4660 0.6740 0.1880],'EdgeColor',[1,1,1])
    l1 = xline(nanmedian(Edist_peA),'r-','linewidth',1.9);
    l0 = xline(nanmean(Edist_peA),'b-','linewidth',1.9);
    l2 = xline(UCI_peA,'k--','linewidth',1.3);
    %xline(quantile(Edist_pe,0.95),'r--','linewidth',1.3);
    %l3 = xline(true_pol,'m--','linewidth',1.9);
    xline(BCI_peA,'k--','linewidth',1.3)
    xlabel('$\gamma$','interpreter','latex','FontSize',opt.fz2)
    ylabel('Frequency','interpreter','latex','FontSize',opt.fz2)
    title(['Nobs: ',num2str(nobs_peA)],'interpreter','latex','FontSize',opt.fz2)
    %xlim([-0.9,0.5])
    %ylim([0,50])
    legend([l0,l1,l2],{'Mean','Median','$90\%$ CI'},'Location','best','FontSize',opt.fz2-5,'interpreter','latex','box','off')
    
    % Around the median Identification window
    figure(24);
    hold on
    set(gca, 'FontSize',opt.fz1);
    histogram(Edist_peB,'BinWidth',0.1/2.5,'FaceColor',[0.4660 0.6740 0.1880],'EdgeColor',[1,1,1])
    l1 = xline(nanmedian(Edist_peB),'r-','linewidth',1.9);
    l0 = xline(nanmean(Edist_peB),'b-','linewidth',1.9);
    l2 = xline(UCI_peB,'k--','linewidth',1.3);
    %xline(quantile(Edist_pe,0.95),'r--','linewidth',1.3);
    %l3 = xline(true_pol,'m--','linewidth',1.9);
    xline(BCI_peB,'k--','linewidth',1.3)
    xlabel('$\gamma$','interpreter','latex','FontSize',opt.fz2)
    ylabel('Frequency','interpreter','latex','FontSize',opt.fz2)
    title(['Nobs: ',num2str(nobs_peB)],'interpreter','latex','FontSize',opt.fz2)
    %xlim([-0.9,0.5])
    %ylim([0,50])
    legend([l0,l1,l2],{'Mean','Median','$90\%$ CI'},'Location','best','FontSize',opt.fz2-5,'interpreter','latex','box','off')
    
    %% Unrestricted 
    
    % Pooled Policy Effect
    figure(25);
    hold on
    set(gca, 'FontSize',opt.fz1);
    histogram(Edist_pe,'BinWidth',0.1/2.5,'FaceColor',[0.4660 0.6740 0.1880],'EdgeColor',[1,1,1])
    l1 = xline(nanmedian(Edist_pe),'r-','linewidth',1.9);
    l0 = xline(nanmean(Edist_pe),'b-','linewidth',1.9);
    l2 = xline(UCI_pe,'k--','linewidth',1.3);
    %xline(quantile(Edist_pe,0.95),'r--','linewidth',1.3);
    %l3 = xline(true_pol,'m--','linewidth',1.9);
    xline(BCI_pe,'k--','linewidth',1.3)
    xlabel('$\gamma$','interpreter','latex','FontSize',opt.fz2)
    ylabel('Frequency','interpreter','latex','FontSize',opt.fz2)
    title(['Nobs: ',num2str(nobs_pe)],'interpreter','latex','FontSize',opt.fz2)
    %xlim([-0.9,0.5])
    %ylim([0,50])
    legend([l0,l1,l2],{'Mean','Median','$90\%$ CI'},'Location','best','FontSize',opt.fz2-5,'interpreter','latex','box','off')
    
    % Identification Interval
    figure(26)
    hold on
    set(gca, 'FontSize',opt.fz1);
    histogram(Edist_tpn-par.tp,'BinWidth',1.1,'FaceColor',[0.4660 0.6740 0.1880],'EdgeColor',[1,1,1])
    l1 = xline(nanmedian(Edist_tpn-par.tp),'r-','linewidth',1.9);
    l0 = xline(nanmean(Edist_tpn-par.tp),'b-','linewidth',1.9);
%    l3 = xline(true_tp-par.tp,'m--','linewidth',1.9);
    l2 = xline(UCI_tpn-par.tp,'k--','linewidth',1.3);
    xline(BCI_tpn-par.tp,'k--','linewidth',1.3)
    xlabel('Size of Identification Window (Stages)','interpreter','latex','FontSize',opt.fz2)
    ylabel('Frequency','interpreter','latex','FontSize',opt.fz2)
    title(['Nobs: ',num2str(nobs_tpn)],'interpreter','latex','FontSize',opt.fz2)
    %xlim([37-par.tp,65-par.tp])
    %ylim([0,280])
    legend([l0,l1,l2],{'Mean','Median','$90\%$ CI'},'Location','best','FontSize',opt.fz2-5,'interpreter','latex','box','off')

    for k=[20:26]   
       saveas(figure(k),sprintf(['output/f%d_',par.name_fapp,par.fapp,'.png'],k)); % will create FIG1, FIG2,...
    end
    end

end
%------------------------------------------------------------------
function flagC = detCT(data,par,opt)
%{ 
This function funtion runs the first lines of norm_func to detect the
leading region
Throws and error if the algorithm cannot detect a leading region
%}

% Pick one to be the control


    par.Cname = char(data.names_st(1,:));
    par.Tname = char(data.names_st(2,:));

    data.C = data.outcome(:,1);
    data.T = data.outcome(:,2);
    
    %% Smoothing
    [datas]=smooth_ts(data,par,opt);
    %% Test Mapping Region 1 to Region 2
    I = datas.C(par.tu)>datas.T(par.tu);
    if I==1
        data.T = datas.T;
        data.C = datas.C;
        data.TO = datas.TO;
        data.CO = datas.CO;
        par.Cname = char(data.names_st(1,:));
        par.Tname = char(data.names_st(2,:));
    else
        data.T = datas.C;
        data.C = datas.T;
        data.TO = datas.CO;
        data.CO = datas.TO;
        par.Cname = char(data.names_st(2,:));
        par.Tname = char(data.names_st(1,:));
    end
      
   
    opt.manual_guess = 2;   % 1: Pick Coefficient Initial Guess Manually 2: Polinomial analitical guess
    opt.sgues = 0;          % Use fmincon with sgues
    par.m_guess = [0.8,10 ,1];
   
    if opt.manual_guess==1
        if size(par.m_guess,2)~=par.nm+1
            error('ERROR: Wrong number of normalization params provided for the guess')
        end
    end
        
    [olp1,flag1,flag_ext1] = norm_func_alt(data,par,opt);
    
 
  
    if flag_ext1==1 
        error('Extreme case detected, try changing the guess')
    end
    
    if flag1<=0 
        error('Minimization didnt converge. Algorithm coudnt detect leading region, debug manually');
    else  
        if olp1>0
            if I ==1
                flagC = 1;
            else
                flagC = 2;
            end
            disp('****************************')
            disp('*** Leading region found ***')
            disp('****************************')
            disp(['*** ',par.Cname,' leads ',par.Tname,' ***'])
            disp('****************************')
        elseif olp1<0
            if I ==1
                flagC = 2;
            else
                flagC = 1;
            end
            disp('****************************')
            disp('*** Leading region found ***')
            disp('****************************')
            disp(['*** ',par.Tname,' leads ',par.Cname,' ***'])
            disp('****************************')
        else
            error('Algorithm coudnt detect leading region, debug manually');
        end
   
    end
    
    
end 
%------------------------------------------------------------------------
function [olp,eflagCT,flag_ext]=norm_func_alt(data,par,opt)
close all

%{
figure(100)
hold on
plot(data.CO,'bx')
plot(data.C,'b--')
plot(data.TO,'rx')
plot(data.T,'r--')
%}

% Keep full series
par.nt_long = length(data.C);
par.time_long = par.time;

% Trim to use prepolicy data for normalization
par.time = par.time(1:par.tu,1);
data.C = data.C(1:par.tu,1);

par.nt = length(data.C);

% Recompute some location parameters according to the trimming
data.T = data.T(1:par.tu);
par.tvec_c      = (1:par.nt)';
par.tvec_c_long = (1:par.nt_long)';



% Initial guess for mapping parameters
%{ 
phi(1): Magnitude Shifter
phi(2): Time Shifter
phi(3): Speed Shifter
%}
%% Mapping C to T
m_guess = par.m_guess;

opt.CT = 1;
[cf.phiCT,x0] = make_x(data,opt,m_guess,par); 
if opt.manual_guess==1
    if size(cf.phiCT,2)~=par.nm+1
        error('ERROR: The number of manual normalization params, is not the correct one')
    end
end
if opt.sgues==1
    [x0,~,~] = fmincon(@(x)func_df(x,data.T,data.C,opt,par),x0,[],[],[],[],[0,-Inf,-Inf],[Inf,Inf,Inf],[],opt.sol_fmin);
end
try
    [xCT,fvalCT,eflagCT] = fminunc(@func_df,x0,opt.sol,data.T,data.C,opt,par);
catch
    xCT = m_guess;
    fvalCT = 0;
    eflagCT = -3; 
end
if eflagCT==0
    disp('WARNING! SOLVER C to T QUIT DUE TO MAX ITER')
end
cf.phiCT = xCT;
%disp(['Coefficients in numerical mapping C to T: ', num2str(cf.phiCT)]);

tu_norm = func_stage(par.tu,cf.phiCT,0,par);
flag_ext = 0;
if tu_norm>2*par.nt
    flag_ext = 1; % Extreme case
end
olp = tu_norm - par.tu;

%func_df(xCT,data.T,data.C,opt,par)

end