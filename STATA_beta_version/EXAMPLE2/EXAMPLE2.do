/*
 --------------------
 EXAMPLE 2:
 --------------------
 Standalone, just run this file.

 This example illustrates the use of the sbinrm.ado  routine described in  
 Aleman,  Busch,  Ludwig,  and  Santaeul`alia-Llopis  (2023) 
 For the complete description of the sbinrm.ado functionalities see the README file.

 -----------------------
 DESCRIPTION: 
 -----------------------

 Exact identification using a polynomial of degree 3 as the data generating
 process.
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
Outcome for region 1 = rest
Outcome for region 2 = reg1
tp = -1  Year of policy implementation
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
import delimited "https://raw.githubusercontent.com/alemanchris/SBI_public/main/STATA_beta_version/EXAMPLE2/input/example2_data.csv" 

* Running sbinrm command
sbinrm  time rest reg1, tp(-1) np(4) sp(0)

* Load saved data
clear all
use dat_NM
qui: replace XN = . if XN==-99999

* Graph Before Normalizatio
graph twoway (line YCS X if X<=-1, sort  lcolor(blue) lwidth(medium)) (line YTS X if X<=-1, sort lcolor(red) lwidth(medium)) (scatter YC X, sort  mcolor(blue) ) (scatter YT X, sort mcolor(red)), xline(-1) legend(on order(3 "yC(t): REST"  4 "yT(t): REG1") size(large) ring(0) bplacement(north)) name(graph_before, replace)  ylabel(0[5] 28, labsize(vlarge)) xtitle("Time", size(vlarge)) xlabel(-5[1]5, labsize(vlarge)) ytitle("My Outcome", size(vlarge))

* Graph After Normalization
graph twoway (line YCNS XN if XN<=-1, sort  lcolor(blue) lwidth(medium)) (line YTS X if X<=-1, sort lcolor(red) lwidth(medium)) (scatter YCN XN if XN<=5, sort  mcolor(blue) ) (scatter YT X, sort mcolor(red)), xline(-1) xline(0.319) legend(on order(3 "yC(s;\psi): REST Norm"  4 "yT(s): REG1") size(large) ring(0) bplacement(north)) name(graph_after, replace)  ylabel(0[5] 28, labsize(vlarge)) xtitle("Stage", size(vlarge)) xlabel(-5[1]5, labsize(vlarge)) ytitle("My Outcome", size(vlarge))

* Graph Policy effect
clear all
use dat_PE
qui: replace X = . if X==-99999
graph twoway (connected Y X if X>=-2, sort  mcolor(black) lcolor(black) lwidth(thick)) , yline(0) xline(-1) xline(0.319) legend(on order(1 "\gamma(s)") size(large) ring(0) bplacement(north)) name(graph_gamma, replace)  ylabel(-0.1[0.05] 0.1, labsize(vlarge)) xtitle("Stage", size(vlarge)) xlabel(-2[1]2, labsize(vlarge)) ytitle("\gamma(s)", size(vlarge))

* Compare analytical results with SBI results
qui: do "https://raw.githubusercontent.com/alemanchris/SBI_public/main/STATA_beta_version/EXAMPLE2/input/get_analytical_sols2"
disp("Compare analytical results with SBI results")
list Coeffs Analytical SBI_norm
