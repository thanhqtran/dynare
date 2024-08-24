    // Model 11. Home production
    // Dynare code
    // File: model10.mod
    // José L. Torres. University of Málaga (Spain)
    // Endogenous variables
    var Y, Cm, Ch, I, K, Lm, Lh, W, R, A, B;
    // Exogenous variables
    varexo e, u;
    // Parameters
    parameters alpha, beta, delta, gamma, omega, eta, theta, rho1, rho2;
    // Calibration
    alpha = 0.35;
    beta  = 0.97;
    delta = 0.06;
    gamma = 0.40;
    omega = 0.45;
    eta   = 0.80;
    theta = 0.80;
    rho1  = 0.95;
    rho2  = 0.95;
    // Equations of the model economy
    model;
    gamma*omega*(Cm^(eta-1))/(omega*Cm^eta+(1-omega)*Ch^eta)
       =(1-gamma)/(W*(1-Lm-Lh));
    gamma*(1-omega)*(Ch^(eta-1))/(omega*Cm^eta+(1-omega)*Ch^eta)
       =(1-gamma)/(B*(1-Lm-Lh));
    ((Cm^(eta-1))/(omega*Cm^eta+(1-omega)*Ch^eta))/((Cm(+1)^(eta-1))
       /(omega*Cm(+1)^eta+(1-omega)*Ch(+1)^eta))=beta*(R(+1)+1-delta);
    Y = A*(K(-1)^alpha)*(Lm^(1-alpha));
    Ch = B*Lh^theta;
    K = (Y-Cm)+(1-delta)*K(-1);
    I = Y-Cm;
    W = (1-alpha)*A*(K(-1)^alpha)*(Lm^(-alpha));
    R = alpha*A*(K(-1)^(alpha-1))*(Lm^(1-alpha));
    log(A) = rho1*log(A(-1))+e;
    log(B) = rho2*log(B(-1))+u;
    end;
    // Initial values
    initval;
    Y = 1; 
    Cm = 0.75;
    Ch = 0.2; 
    Lm = 0.3;
    Lh = 0.1; 
    K = 3.5;
    I = 0.25;
    W = (1-alpha)*Y/Lm;
    R = alpha*Y/K;
    A = 1;
    B = 1;
    e = 0;
    u = 0;
    end;
     // Steady state
    steady;
    // Blanchard-Kahn conditions
    check;
    // Disturbance analysis
    shocks;
    var e; stderr 0.01;
    var u; stderr 0.01;
    end;
     // Stochastic simulation
    stoch_simul;
