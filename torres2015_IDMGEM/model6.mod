// Model 6: Investment Specific Technological Change
// K(t+1)=(1-delta)*K(t)+Z(t)*I(t)
// Dynare code
// File: model6.mod
// José L. Torres. University of Málaga (Spain)
// Endogenous variables
var Y, C, I, K, L, W, R, A, Z;
// Exogenous variables
varexo e, u;
// Parameters
parameters alpha, beta, delta, gamma, rho1, rho2;
// Calibration of the parameters
alpha = 0.35;
beta = 0.97;
delta = 0.06;
gamma = 0.40;
rho1 = 0.95;
rho2 = 0.95;
// Equations of the model economy
model;
C = (gamma/(1-gamma))*(1-L)*(1-alpha)*Y/L;
1 = beta*(Z*C/(Z(+1)*C(+1)))
 *(Z*alpha*Y(+1)/K+(1-delta));
Y = A*(K(-1)^alpha)*(L^(1-alpha));
K = Z*I+(1-delta)*K(-1);
I = Y-C;
W = (1-alpha)*A*(K(-1)^alpha)*(L^(-alpha));
R = alpha*A*(K(-1)^(alpha-1))*(L^(1-alpha));
log(A) = rho1*log(A(-1))+e;
log(Z) = rho2*log(Z(-1))+u;
end;
// Initial values
initval;
Y = 1; 
C = 0.8; 
L = 0.3; 
K = 3.5;
I = 0.2;
W = (1-alpha)*Y/L;
R = alpha*Y/K;
A = 1;
Z = 1;
e = 0;
u = 0;
end;
 // Steady state
steady;
// Blanchard-Kahn conditions
check;
// Perturbation analysis
shocks;
var e; stderr 0.01;
var u; stderr 0.01;
end;
// Stochastic simulation
stoch_simul;
