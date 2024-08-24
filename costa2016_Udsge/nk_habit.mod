// NK model with household rigidity
// (Habit formation and non-Ricardian agents)
// Chapter 5 (UNDERSTANDING DSGE MODELS)

var Y I C CR CNR R K W L LR LNR PIW P PI A LAMBDAR LAMBDANR CM;
varexo e;

parameters sigma phi alpha beta delta rhoa psi theta thetaW psiW phic omegaR;
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

model(linear);
#Pss = 1;
#Rss = Pss*((1/beta)-(1-delta));
#CMss = ((psi-1)/psi)*(1-beta*theta)*Pss;
#Wss = (1-alpha)*(CMss^(1/(1-alpha)))*((alpha/Rss)^(alpha/(1-alpha)));
#Yss = ((Rss/(Rss-delta*alpha*CMss))^(sigma/(sigma+phi))) *((1-phic*beta)*((1-phic)^(-sigma))*(1-beta*thetaW)*((psiW-1)/psiW)*(Wss/Pss) *(Wss/((1-alpha)*CMss))^phi)^(1/(sigma+phi));
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

//1-Ricardian household Lagrangian
LAMBDAR = (sigma/((1-phic)*(1-phic*beta)))*(phic*beta*(CR(+1)-phic*CR) -(CR-phic*CR(-1)))-P;
//2-Phillips equation for Ricardian households’ wages
PIW = beta*PIW(+1)+((1-thetaW)*(1-beta*thetaW)/thetaW)*(LR-LAMBDAR-W); //3-Gross inflation rate for wages
PIW = W - W(-1);
//4-Euler equation
(Pss/beta)*(LAMBDAR+P-LAMBDAR(+1))=(1-delta)*Pss*P(+1) + Rss*R(+1);
//5-Law of motion of capital
K = (1-delta)*K(-1) + delta*I;
//6-Non-Ricardian household Lagrangian
LAMBDANR = (sigma/((1-phic)*(1-phic*beta)))*(phic*beta*(CNR(+1)-phic*CNR) -(CNR-phic*CNR(-1)))-P;
//7-Phillips equation for non-Ricardian households’ wages
PIW = beta*PIW(+1)+((1-thetaW)*(1-beta*thetaW)/thetaW)*(LNR-LAMBDANR-W); 
//8-Non-Ricardian household budget constraint
P+CNR=W+LNR;
//9-Aggregate consumption
Css*C = omegaR*CRss*CR + (1-omegaR)*CNRss*CNR;
//10-Aggregate labor
Lss*L = omegaR*LRss*LR + (1-omegaR)*LNRss*LNR;
//11-Production function
Y = A + alpha*K(-1) + (1-alpha)*L;
//12-Demand for capital
K(-1) = Y - R;
//13-Demand for labor
L = Y - W;
//14-Marginal cost
CM = ((1-alpha)*W + alpha*R - A);
//15-Phillips equation
PI = beta*PI(+1)+((1-theta)*(1-beta*theta)/theta)*(CM-P);
//16-Gross inflation rate
PI = P - P(-1);
//17-Equilibirum condition
Yss*Y = Css*C + Iss*I;
//18-Productivity shock
A = rhoa*A(-1) + e;
end;

model_diagnostics;
steady;
check (qz_zero_threshold=1e-20);
shocks;

var e; 
stderr 0.01; 
end;
stoch_simul(qz_zero_threshold=1e-20) Y I C CR CNR R K W L LR LNR PIW PI A;
