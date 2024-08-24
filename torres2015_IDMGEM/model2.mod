// Model 2: Basic DSGE model
// Dynare code
// File: model2.mod
// José L. Torres. University of Málaga (Spain)
// Endogenous variables
var Y, C, I, K, L, W, R, A;
// Exogenous variables
varexo e;
// Parameters
parameters alpha, beta, delta, gamma, rho;
// Calibration
alpha = 0.35;
beta  = 0.97;
delta = 0.06;
gamma = 0.40;
rho   = 0.95;
// Equations of the model
model;
C = (gamma/(1-gamma))*(1-L)*(1-alpha)*Y/L;
1 = beta*((C/C(+1))*(R(+1)+(1-delta)));
Y = A*(K(-1)^alpha)*(L^(1-alpha));
K = (Y-C)+(1-delta)*K(-1);
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
W = (1-alpha)*Y/L;
R = alpha*Y/K;
A = 1;
e = 0;
end;
// Steady State computation
steady;
// Blanchard-Kahn conditions
check;
// Shock analysis: TFP shock
shocks;
var e; stderr 0.01;
end;
// Stochastic simulation
stoch_simul;
