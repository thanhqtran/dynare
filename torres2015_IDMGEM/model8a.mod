    // Model 8a: Permanent change in public consumption
    // Dynare code
    // File: model8a.mod
    // José L. Torres. University of Málaga (Spain)
    // Fixed figure command issues
    // Endogenous variables
    var Y, Cp, Cg, I, K, L, R, W, A;
    // Exogenous variables
    varexo e, tauc, taul, tauk, zita;
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
    Cg = zita*Y;
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
    zita = 0.100;
    tauc = 0.116;
    tauk = 0.225;
    taul = 0.344;
    end;
     // Steady state
    steady;
    SS0=oo_.steady_state;
    // Blanchard and Khan conditions
    check;
    // Final values
    endval;
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
    zita = 0.120;
    tauc = 0.116;
    tauk = 0.225;
    taul = 0.344;
    end;
     // Steady state
    steady;
    // Perturbation
    shocks;
    var zita;
    // Disturbance periods
    periods 0;
    // Change to final value
    values 0;
    end;
    // Deterministic simulation
    simul(periods=58);
    // Figures
    figure;
    subplot(2,2,1);
    plot(oo_.endo_simul(1,:) - SS0(1));
    title('Output');
    subplot(2,2,2);
    plot(oo_.endo_simul(2,:) - SS0(2));
    title('Private consumption');
    subplot(2,2,3);
    plot(oo_.endo_simul(3,:) - SS0(3));
    title('Public consumption');
    subplot(2,2,4);
    plot(oo_.endo_simul(4,:) - SS0(4));
    title('Investment');
    
    figure;
    subplot(2,2,1);
    plot(oo_.endo_simul(5,:) - SS0(5));  // K should correspond to SS0(5)
    title('Capital stock');
    subplot(2,2,2);
    plot(oo_.endo_simul(6,:) - SS0(6));  // L should correspond to SS0(6)
    title('Working hours');
    subplot(2,2,3);
    plot(oo_.endo_simul(7,:) - SS0(7));  // R should correspond to SS0(7)
    title('Rental rate of capital');
    subplot(2,2,4);
    plot(oo_.endo_simul(8,:) - SS0(8));  // W should correspond to SS0(8)
    title('Wage');
