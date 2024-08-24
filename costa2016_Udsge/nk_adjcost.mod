// NK Model with Rigidity in Investments
// Costs relating to investment adjustment and the 
// non-utilization of maximum installed capacity 
// - Chapter 6 (UNDERSTANDING DSGE MODELS)

var Y I C CR CNR R K W L LR LNR PIW P PI A LAMBDAR LAMBDANR Q U CM; 
varexo e;

parameters sigma phi alpha beta delta rhoa psi theta
thetaW psiW phic omegaR Psi1 Psi2 chi;
sigma = 2;
phi = 1.5; 
alpha = 0.35; 
beta = 0.985; 
delta = 0.025; 
rhoa = 0.95; 
psi = 8;
theta = 0.75;
thetaW = 0.75;
psiW = 21;
phic = 0.8;
omegaR = 0.5;
Psi1 = ((1/beta)-(1-delta)); Psi2 = 1;
chi = 1;

model(linear);

#Uss = 1;
#Pss = 1;
#Rss = Pss*Psi1;
#CMss = ((psi-1)/psi)*(1-beta*theta)*Pss;
#Wss = (1-alpha)*(CMss^(1/(1-alpha)))*((alpha/Rss)^(alpha/(1-alpha))); 
#Yss = ((Rss/(Rss-delta*alpha*CMss))^(sigma/(sigma+phi)))*((1-phic*beta) *((1-phic)^(-sigma))*(1-beta*thetaW)*((psiW-1)/psiW)*(Wss/Pss) *(Wss/((1-alpha)*CMss))^phi)^(1/(sigma+phi));
#Kss = alpha*CMss*(Yss/Rss); 
#Iss = delta*Kss;
#Css = Yss - Iss;
#Lss = (1-alpha)*CMss*(Yss/Wss); 
#CRss = Css;
#CNRss = Css;
#LRss = Lss;
#LNRss = Lss;
#LAMBDARss = (1/Pss)*(CRss^(-sigma))*((1-phic)^(-sigma))*(1-phic*beta); 
#LAMBDANRss = (1/Pss)*(CNRss^(-sigma))*((1-phic)^(-sigma))*(1-phic*beta); 
#Qss = Pss*LAMBDARss;

// 1-Ricardian household Lagrangian
LAMBDAR = (sigma/((1-phic)*(1-phic*beta))) *(phic*beta*(CR(+1)-phic*CR)-(CR-phic*CR(-1)))-P;
// 2-Phillips equation for Ricardian household wages
PIW = beta*PIW(+1)+((1-thetaW)*(1-beta*thetaW)/thetaW)*(LR-LAMBDAR-W); 
// 3-Gross wage inflation rate
PIW = W - W(-1);
// 4-Tobin’s Q
(Qss/beta)*Q = (1-delta)*Qss*Q(+1)+LAMBDARss*Rss*Uss *(LAMBDAR(+1)+R(+1)+U(+1))-LAMBDARss*Pss*Psi1*Uss*U(+1);
// 5-Demand for Installed Capacity
(Rss/Pss)*(R-P) = Psi2*Uss*U;
// 6-Demand for Investment 
LAMBDARss*Pss*(LAMBDAR+P)-Qss*Q+chi*Qss*(I-I(-1))=chi*beta*Qss*(I(+1)-I); 
// 7-Law of Motion of Capital
K = (1-delta)*K(-1) + delta*I;
// 8-Non-Ricardian household Lagrangian
LAMBDANR = (sigma/((1-phic)*(1-phic*beta))) *(phic*beta*(CNR(+1)-phic*CNR)-(CNR-phic*CNR(-1)))-P;
// 9-Phillips equation for non-Ricardian household wages
PIW = beta*PIW(+1)+((1-thetaW)*(1-beta*thetaW)/thetaW)*(LNR-LAMBDANR-W); 
// 10-Budget constraint of the non-Ricardian household
P+CNR=W+LNR;
// 11-Aggregate consumption
Css*C = omegaR*CRss*CR + (1-omegaR)*CNRss*CNR;
// 12-Aggregate labor
Lss*L = omegaR*LRss*LR + (1-omegaR)*LNRss*LNR;
// 13-Production function
Y = A + alpha*(U+K(-1)) + (1-alpha)*L;
// 14-Demand for capital
U + K(-1) = Y - R;
// 15-Demand for labor
L = Y - W;
// 16-Marginal cost
CM = ((1-alpha)*W + alpha*R - A);
// 17-Phillips Equation
PI = beta*PI(+1)+((1-theta)*(1-beta*theta)/theta)*(CM-P);
// 18-Gross inflation rate
PI = P - P(-1);
// 19-Equilibrium condition
Yss*Y = Css*C + Iss*I;
// 20-Productivity shock
A = rhoa*A(-1) + e;
end;

model_diagnostics;
steady;
check (qz_zero_threshold=1e-20);
shocks;

var e; 
stderr 0.01; 
end;

stoch_simul(qz_zero_threshold=1e-20)
Y I C CR CNR R K W L LR LNR PIW PI Q U A;