    // Model 7c. Taxes: Announced transitory change in the consumption tax rate
    // Dynare code
    // File: model7c.mod
    // José L. Torres. University of Málaga (Spain)
    // Endogenous variables
    var Y, C, I, F, K, L, R, W, A;
    // Exogenous variables
    varexo e, tauc, taul, tauk;
    // Parameters
    parameters alpha, beta, delta, gamma, rho;
    // Calibrated parameters
    alpha = 0.35;
    beta  = 0.97;
    delta = 0.06;
    gamma = 0.40;
    rho   = 0.95;
    // Equations of the model
    model;
    (1+tauc)*C=(gamma/(1-gamma))*(1-L)*(1-taul)*(1-alpha)*Y/L;
    1 = beta*((((1+tauc)*C)/((1+tauc(+1))*C(+1)))*((1-tauk)*(R(+1)-delta)+1));
    Y = A*(K(-1)^alpha)*(L^(1-alpha));
    K = I+(1-delta)*K(-1);
    I = Y-C;
    W = (1-alpha)*A*(K(-1)^alpha)*(L^(-alpha));
    R = alpha*A*(K(-1)^(alpha-1))*(L^(1-alpha));
    F = tauc*C+taul*W*L+tauk*(R-delta)*K;
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
    tauc = 0.116;
    tauk = 0.225;
    taul = 0.344;
    F = tauc*C+taul*W*L+tauk*(R-delta)*K;
    end;
     // Steady state
    steady;
    SS0=oo_.steady_state;
    // Blanchard-Kahn conditions
    check;
    // Final values
    endval;
    Y = 1; 
    C = 0.8; 
    L = 0.3; 
    K = 3.5;
    I = 0.2;
    W = (1-alpha)*Y/L;
    R = alpha*Y/K;
    A = 1;
    e = 0;
    tauc = 0.116;
    tauk = 0.225;
    taul = 0.344;
    F = tauc*C+taul*W*L+tauk*(R-delta)*K;
    end;
     // Steady state
    steady;
    // Disturbance: change in the consumption tax
    shocks;
    var tauc;
    // Period of the change
    periods 4:8;
    // Transitory consumption tax rate 
    values 0.13;
    end;
    // Deterministic simulation
    simul(periods=38);
    // Figures
    figure;
    subplot(2,2,1);
    plot(oo_.endo_simul(1,:) - SS0(1));
    title('Output');
    subplot(2,2,2);
    plot(oo_.endo_simul(2,:) - SS0(2));
    title('Consumption');
    subplot(2,2,3);
    plot(oo_.endo_simul(3,:) - SS0(3));
    title('Investment');
    subplot(2,2,4);
    plot(oo_.endo_simul(4,:) - SS0(9));  // F should correspond to SS0(9)
    title('Fiscal revenues');
    
    figure;
    subplot(2,2,1);
    plot(oo_.endo_simul(5,:) - SS0(4));  // K should correspond to SS0(4)
    title('Capital stock');
    subplot(2,2,2);
    plot(oo_.endo_simul(6,:) - SS0(5));  // L should correspond to SS0(5)
    title('Worked hours');
    subplot(2,2,3);
    plot(oo_.endo_simul(7,:) - SS0(6));  // R should correspond to SS0(6)
    title('Interest rate');
    subplot(2,2,4);
    plot(oo_.endo_simul(8,:) - SS0(7));  // W should correspond to SS0(7)
    title('Wage');
