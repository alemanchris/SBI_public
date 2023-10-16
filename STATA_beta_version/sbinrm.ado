program define sbinrm, eclass
*! version 1.0.0 
*! Date November 7th 2022
*! version 2.0.0 
*! Date September 8th 2023
        version 14
 
        syntax varlist, tp(integer) np(integer) sp(integer)
		marksample touse
		fvexpand `varlist' 		
        local cnames `r(varlist)'
        mata tua = `tp'
	mata nup = `np'
	mata smoa = `sp'
	mata SBI_parse("`cnames'", "`touse'", tua, nup, smoa)
end
//------------------------------------------------------------------------------
mata
// define interpolating function 
//-----------------------------------------------------------------------------
function interp1(x,y,xq)
// Careful this is not interpolating zero, fix that!
{
nel1 = rows(xq)
nel2 = rows(x)
yq = J(nel1,1,-99999)

	for (i=1; i<=nel1; i++) {     // for loop 
		diff_aux = x:-xq[i,1]
		ind_sign = sign(diff_aux)
		aux_del = select(ind_sign, ind_sign[.,1]:<0)
		pos = colsum(abs(aux_del))
		ind = 0
    

		  if (pos > 0) { 
				if (pos<nel2){
				ind = 1
				}
				else{			
				}
		   }
		   else{ 
		   }
			
		   if (ind==1){
				y0 = y[pos,1]
				y1 = y[pos+1,1]
				x0 = x[pos,1]
				x1 = x[pos+1,1]
				aux_del2 = y0 + ((y1 - y0)/(x1 - x0)) * (xq[i,1] - x0)
				yq[i,1] =  aux_del2 //&J(1,1,aux_del2)	  
		   }
		   else{
		   		//yq[i,1]=&J(1,1,NULL)
		   }
	} // end for loop
   
return(yq)
}
//-----------------------------------------------------------------------------
function func_stage(time,ppsi)
{
	tl = rows(time)  // lenght
	y = J(tl,1,ppsi[1,1])
	nm = rows(ppsi)+1
    for (i=3; i<=(nm); i++){
        y = y + (time:^(i-2)):*ppsi[i-1,1]
    }
return(y)
}
//-----------------------------------------------------------------------------
function func_dist(todo,pphi,yT,yC,tu,we,fv,g,H){
	np =  cols(pphi)
	if (np==4){
		ppsi2 = pphi[1,3..np]'		
	}
	 else{
		ppsi2 = pphi[1,2..np]'
	}
	
	dtaT = yT
	dtaC = yC
	svec =  range(1,tu,1)
	tvec = func_stage(svec,ppsi2)
	dtaTT = dtaT[1..tu,1]
	dtaCC = dtaC[1..tu,1]
	
	aux = interp1(tvec,dtaCC,svec)
	// Eliminate -99999
	 pred_C_sc = select(aux, aux[.,1]:>-99999)
	 if (np==4){
	 	pred_C_sc = pphi[1,1]:+pphi[1,2]:*pred_C_sc 
	}
	 else{
		
		pred_C_sc = pphi[1,1]:*pred_C_sc 
	}
	
	
	 svec = select(svec, aux[.,1]:>-99999)
	 dtaTT = select(dtaTT, aux[.,1]:>-99999)
	 
	 if (rows(dtaTT)<5) {
		dist = 1e10 
	 }
	 else{
		dist = dtaTT:-pred_C_sc;
	 }

	tl = rows(dist)  
	if (we==1){
		wght = sqrt(range(1,tl,1))		
	}
	 else{
		wght = J(tl,1,1)
	}
	
    wdist = wght:*dist
  
    fv = 0.5:*((wdist')*wdist)

   //return(fv)
}

//-----------------------------------------------------------------------------
function kumsum(x)
{
// take column verctors
	lgth = rows(x)
    y = J(lgth,1,0)
	y[1,1] = x[1,1]
    for (i=2; i<=(lgth); i++){
        y[i,1] = y[i-1,1] + x[i,1]
    }
return(y)
}
//-----------------------------------------------------------------------------
function Tox(n,x)
{
// Manual Chebyshev polinomials
m = rows(x)
T = J(m,1,1)
dT = J(m,1,0)
ddT = J(m,1,0)
for (i=2; i<=(n); i++){

    if (i == 2){
        T = T,x
        dT = dT,J(m,1,1)
        ddT = ddT,J(m,1,0)
    }
    else{
	Taux = 2.0:*x:*T[.,i-1] - T[.,i-2]
        T = T,Taux 
	dTaux = 2.0:*T[.,i-1] + 2.0:*x:*dT[.,i-1] - dT[.,i-2]
        dT = dT,dTaux 
	ddTaux = 2.0:*dT[.,i-1] + 2.0:*dT[.,i-1] + 2.0:*x:*ddT[.,i-1] - ddT[.,i-2]
        ddT = ddT,ddTaux 
    }
}
return(T)
}
//-----------------------------------------------------------------------------
function coll(x_data,y_data,n,cb)
{
// Manual Chebyshev polinomials

np = rows(x_data)
nt = x_data[np]
a = (nt+1)/2
b = a-1
zvec = (x_data:-a):/b
   if (cb == 0){
	    mb = J(nt,1,1)
	    for  (cc=2; cc<=(n); cc++){
		zvec_aux = zvec:^(cc-1)
		mb = mb,zvec_aux
	    }
	    bbeta = invsym((mb'*mb))*(mb'*y_data)  
	    cf= mb*bbeta
    }
    else{
	 cbb = Tox(n,zvec)
	 bbeta = invsym((cbb'*cbb))*(cbb'*y_data)
	 cf= cbb*bbeta
    }
return(cf)
}
//-----------------------------------------------------------------------------
function intycolor(y0,x0,y,x)
{
// Here I integrate using the trapezoidal rule
ys = y0\y
xs = x0\x
nt = rows(xs)
ntr = nt-1  // Number of trapezoids
fx = J(ntr,1,0)

for (i=2; i<=(nt); i++){
        aux_combi = ys[i-1,1]\ys[i,1]
	ysup = max(aux_combi)
	yslo = min(aux_combi)
	// fx[i-1] = (xs[i]-xs[i-1]):*(ys[i-1]):+(xs[i]-xs[i-1]):*(ys[i]-ys[i-1]):/2
	fx[i-1] = (xs[i]-xs[i-1]):*(yslo):+(xs[i]-xs[i-1]):*(ysup-yslo):/2
}
 
return(fx)
}
//-----------------------------------------------------------------------------
function anpol3(tau,time,dataC,dataT)
{
	t0 = J(tau,1,1)
	t1 = time[1..tau,1]
	t2 = time[1..tau,1]:^2
	t3 = time[1..tau,1]:^3
	t4 = time[1..tau,1]:^4
	X = t0,t1,t2,t3
	YC = dataC[1..tau,1]
	YT = dataT[1..tau,1]
	pDC = invsym((X'*X))*(X'*YC)
	pDT = invsym((X'*X))*(X'*YT)
	auxh2 = (pDT[4,1]:/pDC[4,1]):*(pDC[2,1]:-((1:/3):*(pDC[3,1]:^2):/pDC[4,1]))
	auxh4 = (1/3):*((pDT[3]):^2):/pDT[4]:-pDT[2]
	// Get the inverse
	psi1 = (auxh2/(-auxh4))^0.5
	psi0 = (1:/3):*(pDT[3]:/pDT[4]:*psi1):-(1:/3):*(pDC[3]:/pDC[4])
	alt_w1 = pDT[4]:/(pDC[4]:*psi1:^3)
	alt_w0 = pDT[1]-(alt_w1*(pDC[1]+pDC[2]*psi0+pDC[3]*psi0^2+pDC[4]*psi0^3))
	om0 = alt_w0
	om1 = alt_w1
	psi0 = (-psi0/psi1)
	psi1 = 1/psi1
	phi = J(4,1,1)
	phi[1,1] = om0
	phi[2,1] = om1
	phi[3,1] = psi0
	phi[4,1] = psi1
return(phi)
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//**************************************************************
// CORE ROUTINE STARTS HERE
//**************************************************************
void SBI_parse(string scalar indepvars,  string scalar touse, tua, nup, smoa){
		// Fetch data
		real matrix data_mat
		data_mat  = st_data(., indepvars, touse)

		// Smoother indicator	1: Smooth 0:No Smoothing			
		smo =  smoa

		// Get time
		time_o = data_mat[.,1]
                nt = rows(time_o)
		time = range(1,nt,1)


		// Adjust time 	
		
		//min_time = time_o[1,1]
		//time_adjust = (min_time)-1
		//time = time_o:-time_adjust


		// Adjust Policy date
		//tu =  tua:-time_adjust

		I = time_o:<=tua
		tu = sum(I)   
		
		// Rename outcome	
		yC = data_mat[.,2]
		yT = data_mat[.,3]

		yC_aux = yC[1..tu,1]
		yT_aux = yT[1..tu,1]
		time_aux = time[1..tu,1]


		if (smo==1){
			yC_smo = coll(time_aux,yC_aux,4,1)
			yT_smo = coll(time_aux,yT_aux,4,1)
		}
		else{
			yC_smo = yC_aux
			yT_smo = yT_aux
		}

	
		tup1 = tu+1



		yC_aux1 = yC[tup1..nt,1]
		yT_aux1 = yT[tup1..nt,1]

	

		yC_smo = yC_smo\yC_aux1
		yT_smo = yT_smo\yT_aux1


		



		
		nm = nup  //2						// choose number of param 2: linear 3 quadratic
		we = 0			// Weighted minimization



		
		// Initial guess 
 		if (nm==4){
			x0 = (2.78,0.27,13.38,0.914)
			x0 = (2.5,0.1,10,0.9)
 			x0 = anpol3(tu,time,yC_smo,yT_smo)

			if (x0[1,1]==.){

				x0 = J(1,4,1)
				x0[1,1] = 0
				aux1 = max(yT_smo[.,1])
				aux2 = max(yC_smo[.,1])

				x0[1,2] = aux1/aux2
				I = yT_smo:<=yC_smo[tu,1]
				pos_del = sum(I) 
				aux3 = tu-pos_del
				x0[1,3] = abs(aux3)*(3/4)
				x0[1,4] = 1
			}
			else{
				x0 = x0'
			}

			if (x0[1,2]<=0){
				x0[1,1] = 0
				aux1 = max(yT_smo[.,1])
				aux2 = max(yC_smo[.,1])
				x0[1,2] = aux1/aux2
				I = yT_smo:<=yC_smo[tu,1]
				pos_del = sum(I) 
				aux3 = tu-pos_del
				x0[1,3] = abs(aux3)*(3/4)
				x0[1,4] = 1

			}
			else{
				x0 = x0
			}
		
		}
		 else{
			x0 = (0.9,0,1.1)
		}

		
		 // Perform minimization
		 s=optimize_init()
		 optimize_init_which(s, "min")
		 optimize_init_evaluator(s,&func_dist())
		 optimize_init_argument(s, 1, yT_smo)
		 optimize_init_argument(s, 2, yC_smo)
		 optimize_init_argument(s, 3, tu)
		 optimize_init_argument(s, 4, we)
		 optimize_init_conv_ptol(s, 1e-16)
		 optimize_init_conv_vtol(s, 1e-16)
		 optimize_init_conv_maxiter(s, 2000)
		 optimize_init_params(s,x0)
		 cpphi=optimize(s)
		 // 
		

	
		 if (nm==4){
			cppsi = cpphi[1,3..nm]'		
		}
		 else{
			cppsi = cpphi[1,2..nm]'
		}
		tu_norm = func_stage(tu,cppsi)	   		 
		svec = time
		tvec = func_stage(svec,cppsi)

	
		
		// time_long_y = interp1(par.tvec_c_long,par.time_long,tvec,'linear',NaN); // Convert to time units
		if (nm==4){
			
			CO_norm = cpphi[1,1]:+cpphi[1,2]:*interp1(tvec,yC,svec)
			CO_norm_no = cpphi[1,1]:+cpphi[1,2]:*yC
			COS_norm_no = cpphi[1,1]:+cpphi[1,2]:*yC_smo
		}
		else{			
			CO_norm = cpphi[1,1]:*interp1(tvec,yC,svec)
			CO_norm_no = cpphi[1,1]:*yC
			COS_norm_no = cpphi[1,1]:*yC_smo
		}
		TO_int = interp1(svec,yT,tu_norm)
		
		// Eliminate the -99999 for the computation of gamma, but that messes up the order, figure it out
		// I am not taking into account the half days, do it!
		// Trim on indetification window.
		tu1 = floor(tu_norm)
		C_del =  CO_norm[tu..tu1,1]
		T_del =  yT[tu..tu1,1]
		time_del =  svec[tu..tu1,1]
		aux_del = rows(T_del)
		

	
		unb=intycolor(C_del,time_del,CO_norm_no[tu],tu_norm)
		unr=intycolor(T_del,time_del,TO_int,tu_norm)



		ng = rows(unb)


		unr = kumsum(unr) // rougth approximation for area under red
		unb = kumsum(unb) // rougth approximation for area under blue

		
		gamma_s = J(ng,1,0)
		time_s  = J(ng,1,0)


	
		
		for (i=1; i<=(ng); i++){			
			gamma_s[i,1] =  (unr[i,1]-unb[i,1]):/unb[i,1]
			
			if (i<ng){
			time_s[i,1]  =  time_del[i+1,1]
			}
			else{ 
			time_s[i,1]  =  tu_norm
			}
		}		

	
		

		// Policy Effect
		gamma_val = gamma_s[ng]
		
		// Identification window size
		olp = tu_norm-tu
		
		// Print output

		if (nm==4){
			
			om0 = strofreal(cpphi[1,1])
			om1 = strofreal(cpphi[1,2])
			psi0 = strofreal(cpphi[1,3])
			psi1 = strofreal(cpphi[1,4])
		}
		else{	
			om0 = strofreal(0)	
			om1 = strofreal(cpphi[1,1])
			psi0 = strofreal(cpphi[1,2])
			psi1 = strofreal(cpphi[1,3])
		}


		tu_norm_adjust = interp1(time,time_o,tu_norm)

		gamma_val = strofreal(gamma_val[1,1])
		olp_s = strofreal(olp[1,1])
		stp   = strofreal(tu_norm_adjust)


		s00 = ("omega_0 =",om0)		
		s0 = ("omema_1 =",om1)
		s1 = ("psi_0 =",psi0)
		s2 = ("psi_1 =",psi1)
		s00 = invtokens(s00, " ")
		s0 = invtokens(s0, " ")
		s1 = invtokens(s1, " ")
		s2 = invtokens(s2, " ")
		
		s3 = ("gamma =",gamma_val)
		s4 = ("olp =",olp_s)
		s5 = ("s(tp) =",stp)
		s3 = invtokens(s3, " ")
		s4 = invtokens(s4, " ")
		s5 = invtokens(s5, " ")
		
		s0L = ("********")
		s0LL = ("******************************************")
		s1L = ("Results:")
		s2L = ("Normalization Coefficients:")
		printf("\n%s\n",s0L)
		printf("%s\n",s1L)
		printf("%s\n",s0L)
		printf("\n%s\n",s2L)
		printf("%s\n",s00)
		printf("%s\n",s0)
		printf("%s\n",s1)
		printf("%s\n",s2)
		s3L = ("Policy Effect:")
		printf("\n%s\n",s3L)
		printf("%s\n",s3)
		s4L = ("Length of identification window:")
		printf("\n%s\n",s4L)
		printf("%s\n",s4)
		s5L = ("Position of normalized policy date:")
		printf("\n%s\n",s5L)
		printf("%s\n",s5)
		printf("%s\n","")

	// Save individual data sets to generate graphs in stata
	// Before and after Normalization

	ldta = rows(time)
	time_nn = time
	time = time_o
	tvec = interp1(time_nn,time_o,tvec)
	//time = time:+time_adjust
	//tvec = tvec:+time_adjust
	stata("clear")
	st_addobs(ldta)
	(void) st_addvar("float","YCS")
	st_store(.,"YCS",yC_smo)
	(void) st_addvar("float","YC")
	st_store(.,"YC",yC)
	(void) st_addvar("float","YCNS")
	st_store(.,"YCNS",COS_norm_no)
	(void) st_addvar("float","YCN")
	st_store(.,"YCN",CO_norm_no)
	(void) st_addvar("float","YTS")
	st_store(.,"YTS",yT_smo)
	(void) st_addvar("float","YT")
	st_store(.,"YT",yT)
	(void) st_addvar("float","X")
	st_store(.,"X",time)
	(void) st_addvar("float","XN")
	st_store(.,"XN",tvec)
	stata("save dat_NM, replace")
	stata("drop X")
	stata("drop XN")
	stata("drop YCS")
	stata("drop YC")
	stata("drop YCNS")
	stata("drop YCN")
	stata("drop YT")
	stata("drop YTS")


	// Cummulative Policy Effect \gamma(s)

	ldta = rows(gamma_s)
	//time_s = time_s:+time_adjust
	time_s = interp1(time_nn,time_o,time_s)
	stata("clear")
	st_addobs(ldta)
	(void) st_addvar("float","Y")
	st_store(.,"Y",gamma_s)
	(void) st_addvar("float","X")
	st_store(.,"X",time_s)
	stata("save dat_PE, replace")
	stata("drop X")
	stata("drop Y")

	
	



}
end

