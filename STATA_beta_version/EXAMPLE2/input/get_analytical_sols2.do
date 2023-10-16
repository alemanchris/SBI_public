/*
Exact identification using a polynomial of degree 3 as the data generating process.
*/

mata

// Polynomial values
pDC= J(1,4,1)
pDC= 16.019921874999987, 1.023046875000003, -0.043359375000000, 0.000390625000000
pDT= J(1,4,1)
pDT= 0.148901367187497, 0.775327148437501, -0.020170898437500,  0.000142382812500 

// Analytical solution to the normalization  
  
auxh2 = (pDT[4]/pDC[4])*(pDC[2]-((1/3)*(pDC[3]^2)/pDC[4]))
auxh4 = (1/3)*(pDT[3])^2/pDT[4]-pDT[2]

psi1 = (auxh2/(-auxh4))^0.5
psi0 = (1/3)*(pDT[3]/pDT[4]*psi1)-(1/3)*(pDC[3]/pDC[4])
alt_w1 = pDT[4]/(pDC[4]*psi1^3)
alt_w0 = pDT[1]-(alt_w1*(pDC[1]+pDC[2]*psi0+pDC[3]*psi0^2+pDC[4]*psi0^3))

// Get the inverse
phi= J(4,1,1)
phi[1,1] = alt_w0  // om0
phi[2,1] = alt_w1  // om1
phi[3,1] = (-psi0/psi1)  // psi0
phi[4,1] = 1/psi1         // psi1


ldta = rows(phi)
stata("clear")
st_addobs(ldta)
(void) st_addvar("float","Analytical")
st_store(.,"Analytical",phi)
	
end
set obs 4
gen Coeffs = ""
replace Coeffs = "omega_0" in 1
replace Coeffs = "omega_1" in 2
replace Coeffs = "psi_0" in 3
replace Coeffs = "psi_1" in 4
gen SBI_norm = .
replace SBI_norm = 2.792089 in 1
replace SBI_norm =  0.2786751 in 2
replace SBI_norm =  13.38801 in 3
replace SBI_norm =  0.9141033 in 4



