    // Model 8b: Stochastic change in public consumption
    // Dynare code
    // File: model8b.mod
    // José L. Torres. University of Málaga (Spain)
    // remove redundant shock zita
    // Endogenous variables
    var Y, Cp, Cg, I, K, L, R, W, A;
    // Exogenous variables
    varexo e, u, tauc, taul, tauk;
    // Parameters
    parameters alpha, beta, delta, gamma, pi, rho;
    // Calibration of parameters
    alpha = 0.35;
    beta  = 0.97;
    delta = 0.06;
    gamma = 0.40;
    pi    = 0.50;
    rho   = 0.95;
    // Equations of the model economy
    model;
    (1+tauc)*(Cp+pi*Cg)=(gamma/(1-gamma))*(1-L)*(1-taul)*W;
    1 = beta*((1+tauc)*(Cp+pi*Cg)/((1+tauc)*(Cp(+1)+pi*Cg(+1)))
    *((1-tauk)*(R(+1)-delta)+1));
    Y = A*(K(-1)^alpha)*(L^(1-alpha));
    K = I+(1-delta)*K(-1);
    I = Y-Cp-Cg;
    W = (1-alpha)*A*(K(-1)^alpha)*(L^(-alpha));
    R = alpha*A*(K(-1)^(alpha-1))*(L^(1-alpha));
    Cg = u*Y;
    log(A) = rho*log(A(-1))+ e;
    end;
    // Initial values
    initval;
    Y = 1; 
    Cp = 0.8;
    Cg = 0.1; 
    L = 0.3; 
    K = 3.5;
    I = 0.2;
    W = (1-alpha)*Y/L;
    R = alpha*Y/K;
    A = 1;
    e = 0;
    u = 0.1;
    tauc  = 0.116;
    tauk  = 0.225;
    taul  = 0.344;
    end;
    // Steady State
    steady;
    // Blanchard-Kahn condition
    check;
    // Disturbance analysis
    shocks;
    var u; stderr 0.01;
    end;
    // Stochastic simulation
    stoch_simul;
