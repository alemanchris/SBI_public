/*
 --------------------
 EXAMPLE 3:
 --------------------
 Standalone, just run this file.

 This example illustrates the use of the sbinrm.ado  routine described in  
 Aleman,  Busch,  Ludwig,  and  Santaeul`alia-Llopis  (2023) 
 For the complete description of the sbinrm.ado functionalities see the README file.

 -----------------------
 DESCRIPTION: 
 -----------------------

 Exact identification using GLF(Generalized Logistic Func) with as the data generatig process
 This also serves as a placebo example, the policy effect is zero thus
 estimated effect should be zero.
 
 ---------------------------
SYNTAX of sbinrm 
---------------------------
sbirnrm varlist, options
 
varlist: 
- time: time variable (in years, months, etc)
- out1: Outcome for region 1
- out2: Outcome for region 2 

options: 
- tp: time of policy implementation 
- np: number of normalization parameters, either 3 or 4
- sp: smoothing indicator, 1: Smooth (pre-policy) outcome or 0: No smoothing

In this example: 
Outcome for region 1 = BLUE
Outcome for region 2 = RED
tp = 70  Year of policy implementation
np = 4  linear time stage transform (\psi_0, \psi_1) and aditive and proportional level shifter (\omega_0, \omega_1)
sp = 0  Smooth pre-policy outcome 1: Smooth 0: Dont Smooth


*/

* Housekeeping
clear all
set more off
clear programs

* Installing sbinrm routine
qui:net from "https://raw.githubusercontent.com/alemanchris/SBI_public/main/STATA_beta_version/"
net install sbinrm, replace
qui:net from "http://www.stata.com/"

* Fetch Data
import delimited "https://raw.githubusercontent.com/alemanchris/SBI_public/main/STATA_beta_version/EXAMPLE3/input/example3_data.csv" 

* Running sbinrm command
sbinrm  time blue red, tp(70) np(4) sp(0)

* Load saved data
clear all
use dat_NM
qui: replace XN = . if XN==-99999

* Graph Before Normalization
graph twoway (line YCS X if X<=70, sort  lcolor(blue) lwidth(medium)) (line YTS X if X<=70, sort lcolor(red) lwidth(medium)) (scatter YC X, sort  mcolor(blue) ) (scatter YT X, sort mcolor(red)), xline(70) legend(on order(3 "yC(t): BLUE"  4 "yT(t): RED") size(large) ring(0) bplacement(north)) name(graph_before, replace)  ylabel(0[10] 45, labsize(vlarge)) xtitle("Time", size(vlarge)) xlabel(0[30]210, labsize(vlarge)) ytitle("My Outcome", size(vlarge))

* Graph After Normalization
graph twoway (line YCNS XN if XN<=70, sort  lcolor(blue) lwidth(medium)) (line YTS X if X<=70, sort lcolor(red) lwidth(medium)) (scatter YCN XN if XN<=210, sort  mcolor(blue) ) (scatter YT X, sort mcolor(red)), xline(70) xline(82.51) legend(on order(3 "yC(s;\psi): BLUE Norm"  4 "yT(s): RED") size(large) ring(0) bplacement(north)) name(graph_after, replace)  ylabel(0[10] 45, labsize(vlarge)) xtitle("Stage", size(vlarge)) xlabel(0[30]210, labsize(vlarge)) ytitle("My Outcome", size(vlarge))

* Graph Policy effect
clear all
use dat_PE
qui: replace X = . if X==-99999
graph twoway (connected Y X if X>=65, sort  mcolor(black) lcolor(black) lwidth(thick)) , yline(0) xline(70) xline(82.51) legend(on order(1 "\gamma(s)") size(large) ring(0) bplacement(north)) name(graph_gamma, replace)  ylabel(-0.1[0.05] 0.1, labsize(vlarge)) xtitle("Stage", size(vlarge)) xlabel(65[5]90, labsize(vlarge)) ytitle("\gamma(s)", size(vlarge))

* Compare analytical results with SBI results
qui: do "https://raw.githubusercontent.com/alemanchris/SBI_public/main/STATA_beta_version/EXAMPLE3/input/get_analytical_sols3"
disp("Compare analytical results with SBI results")
list Coeffs Analytical SBI_norm



