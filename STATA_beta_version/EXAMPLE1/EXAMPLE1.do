/*
 --------------------
 EXAMPLE 1:
 --------------------
 Standalone, just run this file.

 This example illustrates the use of the sbinrm.ado  routine described in  
 Aleman,  Busch,  Ludwig,  and  Santaeul`alia-Llopis  (2023) 
 For the complete description of the sbinrm.ado functionalities see the README file.

 -----------------------
 DESCRIPTION: 
 -----------------------
 Uses example data to run SBI


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
Outcome for region 1 = reg1
Outcome for region 2 = rest
tp = 1960  Year of policy implementation
np = 3  linear time stage transform (\psi_0, \psi_1) and proportional level shifter (\omega_1)
sp = 1  Smooth pre-policy outcome 1: Smooth 0: Dont Smooth

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
import delimited "https://raw.githubusercontent.com/alemanchris/SBI_public/main/STATA_beta_version/EXAMPLE1/input/example1_data.csv" 

* Running sbinrm command
sbinrm  time reg1 rest, tp(1960) np(3) sp(1)

* Load saved data
clear all
use dat_NM
qui: replace XN = . if XN==-99999

* Graph Before Normalization
graph twoway (line YCS X if X<=1960, sort  lcolor(blue) lwidth(medium)) (line YTS X if X<=1960, sort lcolor(red) lwidth(medium)) (scatter YC X, sort  mcolor(blue) ) (scatter YT X, sort mcolor(red)), xline(1960) legend(on order(3 "yC(t): REG1"  4 "yT(t): REST") size(large) ring(0) bplacement(north)) name(graph_before, replace)  ylabel(0.05[0.04] 0.3, labsize(vlarge)) xtitle("Year", size(vlarge)) xlabel(1934[10]1974, labsize(vlarge)) ytitle("My Outcome", size(vlarge))

* Graph After Normalization
graph twoway (line YCNS XN if XN<=1960, sort  lcolor(blue) lwidth(medium)) (line YTS X if X<=1960, sort lcolor(red) lwidth(medium)) (scatter YCN XN if XN<=1974, sort  mcolor(blue) ) (scatter YT X, sort mcolor(red)), xline(1960) xline(1970.502) legend(on order(3 "yC(s;\psi): REG1 Norm"  4 "yT(s): REST") size(large) ring(0) bplacement(north)) name(graph_after, replace)  ylabel(0.05[0.04] 0.3, labsize(vlarge)) xtitle("Stage", size(vlarge)) xlabel(1934[10]1974, labsize(vlarge)) ytitle("My Outcome", size(vlarge))

* Graph Policy effect
clear all
use dat_PE
qui: replace X = . if X==-99999
graph twoway (connected Y X if X>=1954, sort  mcolor(black) lcolor(black) lwidth(thick)) , yline(0) xline(1960) xline(1970.502) legend(on order(1 "\gamma(s)") size(large) ring(0) bplacement(north)) name(graph_gamma, replace)  ylabel(0[0.1] 0.25, labsize(vlarge)) xtitle("Stage", size(vlarge)) xlabel(1954[5]1974, labsize(vlarge)) ytitle("\gamma(s)", size(vlarge))





