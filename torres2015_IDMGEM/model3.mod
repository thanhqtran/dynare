// Model 3: Habits formation: U[C(t)-phi*C(t-1),O(t)]
// Dynare code
// File: model3.mod
// José L. Torres. University of Málaga (Spain)
// Endogenous variables
var Y, C, I, K, L, W, R, A;
// Exogenous variables
varexo e;
// Parameters
parameters alpha, beta, delta, gamma, rho, phi;
// Calibration
alpha = 0.35;
beta  = 0.97;
delta = 0.06;
gamma = 0.40;
rho   = 0.95;
phi   = 0.80;
// Equations of the model
model;
(gamma/(C-phi*C(-1))-beta*gamma*phi/(C(+1)-phi*C))
    =(1-gamma)/((1-L)*(1-alpha)*Y/L);
(gamma/(C-phi*C(-1))-beta*gamma*phi/(C(+1)-phi*C))/
   (gamma/(C(+1)-phi*C)-beta*gamma*phi/(C(+2)-phi*C(+1)))
    =beta*(alpha*Y(+1)/K+(1-delta));
Y = A*(K(-1)^alpha)*(L^(1-alpha));
K = I+(1-delta)*K(-1);
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
// Steady State
steady;
// Blanchard-Kahn conditions
check;
// Perturbation analysis
shocks;
var e; stderr 0.01;
end;
// Stochastic simulation
stoch_simul;
