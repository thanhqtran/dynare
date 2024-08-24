// Model 5: Investment adjustment costs
// Dynare code
// File: model5.mod
// José L. Torres. University of Málaga (Spain)
// Endogenous variables
var Y, C, I, K, L, W, R, q, A;
// Exogenous variables
varexo e;
// Parameters
parameters alpha, beta, delta, gamma, psi, rho;
// Calibration of the parameters
alpha = 0.35;
beta  = 0.97;
delta = 0.06;
gamma = 0.40;
psi   = 2.00;
rho   = 0.95;
// Equations of the model economy
model;
C=(gamma/(1-gamma))*(1-L)*W;
q=beta*(C/C(+1))*(q(+1)*(1-delta)+R(+1));
q-q*psi/2*((I/I(-1))-1)^2-q*psi*((I/I(-1))-1)*I/I(-1)
 +beta*C/C(+1)*q(+1)*psi*((I(+1)/I)-1)*(I(+1)/I)^2=1;
Y = A*(K(-1)^alpha)*(L^(1-alpha));
K = (1-delta)*K(-1)+(1-(psi/2*(I/I(-1)-1)^2))*I;
I = Y-C;
W = (1-alpha)*A*(K(-1)^alpha)*(L^(-alpha));
R = alpha*A*(K(-1)^(alpha-1))*(L^(1-alpha));
log(A) = rho*log(A(-1))+ e;
end;
// Initial values
initval;
Y = 1; 
C = 0.8; 
L = 0.3; 
K = 3.5;
I = 0.2;
q = 1;
W = (1-alpha)*Y/L;
R = alpha*Y/K;
A = 1;
e = 0;
end;
 // Steady state
steady;
// Blanchard-Kahn conditions
check;
// Perturbation analysis
shocks;
var e; stderr 0.01;
end;
// Stochastic simulation
stoch_simul;
