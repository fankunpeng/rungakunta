function b = laser(x, b, theta, eDelayIndex)
    i = 0;
    a = [];
    t = [];
    h = 5.0e-12;
    transient = 5000.0e-9; 
    tMax = 20.0e-9; 
    div = 10; 
    trans = floor(transient / h);
    n = floor(tMax / h);
    a = [1.3e10, 0.0, 1.90e24];
    r3 = 0.01; % reflectivity of external mirror
    JRatio = 1.11;
    % normalized injection current
    L = 0.225;
    % external cavity length (one-way) [m]
    % fixed parameter values
    r2 = 0.556;
    % reflectivity of internal mirror
    tauIn = 8.0e-12;
    % internal cavity round-trip time [s]
    lambda = 1.537e-6;
    % optical wavelength [m]

    Gn = 8.4e-13; % gain coefficient
    N0 = 1.4e24;
    % carrier density at transparency
    tauP = 1.927e-12; % photon lifetime [s]
    tauS = 2.04e-9; % carrier lifetime [s]
    alpha = 3.0; % alpha parameter
    Nth = 0.0
    Jth = 0.0;
    J = 0.0;
    kap = 0.0;
    tau = 0.0 ;
    detun = 0.0;
    omega0 = 0.0;
    phaseShift = 0.0;
    M = 3;
    C = 2.99792458e8;
    Nth = N0 + 1.0 / tauP / Gn;
    % carrier density at threshold
    Jth = Nth / tauS;
    % injection current at threshold
    J = JRatio * Jth; % injection current
    kap = (1 - r2 * r2) * r3 / r2 / tauIn;
    % injection strength
    tau = 2.0 * L / C;

    b(1) = 1.0 / 2.0 * (Gn * (x(3) - N0) - 1.0 / tauP) * x(1) + kap * eDelayIndex * cos(theta);
    b(2) = alpha / 2.0 * (Gn * (x(2) - N0) - 1.0 / tauP) - kap * eDelayIndex / x(1) * sin(theta);
    b(3) = J - x(3) / tauS - Gn * (x(3) - N0) * x(1) * x(1);
end
