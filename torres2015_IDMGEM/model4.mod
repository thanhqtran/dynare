    // Model 4: Ricardian and non-Ricardian agents
    // Dynare code
    // File: model4.mod
    // José L. Torres. University of Málaga (Spain)
    // Endogenous variables
    var Y, C, C1, C2, I, I1, K, K1, L, L1, L2, W, R, A;
    // Exogenous variables
    varexo e;
    // Parameters
    parameters alpha, beta, delta, gamma, omega, rho;
    // Calibration of the parameters
    alpha = 0.35;
    beta  = 0.97;
    delta = 0.06;
    gamma = 0.40;
    omega = 0.50;
    rho   = 0.95;
    // Equations of the model economy
    model;
    C1=(gamma/(1-gamma))*(1-L1)*W;
    C2=(gamma/(1-gamma))*(1-L2)*W;
    C2=W*L2;
    C =omega*C1+(1-omega)*C2; 
    1 = beta*((C1/C1(+1))
     *(R(+1)+(1-delta)));
    K = omega*K1;
    L = omega*L1+(1-omega)*L2;
    Y = A*(K(-1)^alpha)*(L^(1-alpha));
    K1= I1+(1-delta)*K1(-1);
    I1= W*L1+R*K1-C1;
    I = omega*I1;
    W = (1-alpha)*A*(K(-1)^alpha)*(L^(-alpha));
    R = alpha*A*(K(-1)^(alpha-1))*(L^(1-alpha));
    log(A) = rho*log(A(-1))+ e;
    end;
    // Initial values
    initval;
    Y = 1; 
    C = 0.8;
    C1= 0.6;
    C2= 0.2; 
    L = 0.3;
    L1= 0.3;
    L2= 0.3; 
    K = 3.5;
    K1= 4;
    I = 0.2;
    I1= 0.3;
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
