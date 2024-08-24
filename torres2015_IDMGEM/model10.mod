// Model 10. Human capital
// Dynare code
// File: model10.mod
// José L. Torres. University of Málaga (Spain)
// Endogenous variables
var Y, C, IK, K, IH, H, L, E, W, R, A, B;
// Exogenous variables
varexo e, v;
// Parameters
parameters alpha, beta, deltak, deltah, gamma, theta, rhoA, rhoB;
// Calibration
alpha  = 0.35;
beta   = 0.97;
deltak = 0.06;
deltah = 0.01;
gamma  = 0.40;
theta  = 0.80;
rhoA   = 0.95;
rhoB   = 0.95;
// Equations of the model
model;
C = (gamma/(1-gamma))*(1-L-E)*H*W;
1 = beta*((C/C(+1))*(R(+1)+(1-deltak)));
Y = A*(K(-1)^alpha)*((L*H)^(1-alpha));
K = (Y-C)+(1-deltak)*K(-1);
IK = Y-C;
H=IH+(1-deltah)*H(-1);
IH = B*(E)^theta;
(1-gamma)/((1-L-E)*theta*B*(E)^(theta-1))=
beta*((gamma*W(+1)*L(+1))/C(+1)+((1-gamma)*(1-deltah)/(1-L(+1)-E(+1)*theta*B*(E+1)^(theta-1))));
W = (1-alpha)*A*(K(-1)^alpha)*((L*H)^(-alpha));
R = alpha*A*(K(-1)^(alpha-1))*((L*H)^(1-alpha));
log(A) = rhoA*log(A(-1))+ e;
log(B) = rhoB*log(B(-1))+ v;
end;
// Initial values
initval;
Y  = 1; 
C  = 0.8; 
L  = 0.3; 
K  = 3.5;
IK = 0.2;
K  = 3.5;
E  = 0.15;
IK = 0.15^0.8;
H  = IK/deltah;
W  = (1-alpha)*Y/L;
R  = alpha*Y/K;
A  = 1;
B  = 1;
e  = 0;
v  = 0;
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
