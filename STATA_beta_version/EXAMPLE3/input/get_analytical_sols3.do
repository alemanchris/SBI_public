/*
Exact identification using a polynomial of degree 3 as the data generating process.
*/

mata
// GLF Values for C
theta0C = 0
theta1C = 40
theta2C = 6.93
theta3C = 0.09

// GLF Values for T
theta0T = 0
theta1T = 35
theta2T = 7.92
theta3T = 0.088
// Analytical solution to the normalization  
  
an_psi2 = (theta0T-theta1T)/(theta0C-theta1C);
an_psi1 = theta0T-theta0C*an_psi2
an_psi3 = (theta2C-theta2T)/theta3C;
an_psi4 = theta3T/theta3C;


// Get the inverse
phi= J(4,1,1)
phi[1,1] = an_psi1  // om0
phi[2,1] = an_psi2  // om1
phi[3,1] = (-an_psi3/an_psi4) // psi0
phi[4,1] = 1/an_psi4        // psi1


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
replace SBI_norm = 0.0005213 in 1
replace SBI_norm =  0.8545349 in 2
replace SBI_norm =  11.0661 in 3
replace SBI_norm =  1.020628 in 4



