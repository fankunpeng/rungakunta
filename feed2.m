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
M = 3;
C = 2.99792458e8;

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
plotCnt = 50;
plotMax = 50;
paraMin = 0.0;
paraMax = 0.06;
paraNum = 200;
paraDiv = (paraMax - paraMin) / paraNum;
prev = 0;
now = 0;
next = 0;
Nth = N0 + 1.0 / tauP / Gn;
% carrier density at threshold
Jth = Nth / tauS;
% injection current at threshold
J = JRatio * Jth; % injection current
% injection strength
tau = 2.0 * L / C;
% external-cavity round-trip time for drive
omega0 = 2.0 * pi * C / lambda;
% optical angular frequency
phaseShift = rem(omega0 * tau, 2.0 * pi);
% initial phase shift
delayNum = floor(tau / h);
delayIndex = 1;
DELAY_MAX = 100000;
eDelay = zeros(1, DELAY_MAX);
phiDelay = zeros(1, DELAY_MAX);
arrI = [];
arrA = [];

for i = 1 : DELAY_MAX
    eDelay(i) = a(1);
    phiDelay(i) = a(2);
end

for p = 1 : paraNum 
    r3 = paraMin + paraDiv * p; % reflectivity of external mirror
    kap = (1 - r2 * r2) * r3 / r2 / tauIn;
    
    for i = 1 : trans
        t = h * i;
        % rungakota
        x = zeros(1, M);
        b = zeros(4, M);
        theta = rem(phaseShift + a(2) - phiDelay(delayIndex), 2.0 * pi);
        for i = 1 : 4
            for j = 1 : M
                if i == 1
                    x(j) = a(j);
                elseif i == 2
                    x(j) = a(j) + h * b(1,j) / 2.0;
                elseif i == 3
                    x(j) = a(j) + h * b(2,j) / 2.0;
                elseif i == 4
                    x(j) = a(j) + h * b(3,j) / 2.0;
                end
                b(i,1) = 1.0 / 2.0 * (Gn * (x(3) - N0) - 1.0 / tauP) * x(1)+ kap * eDelay(delayIndex) * cos(theta);
                b(i,2) = alpha / 2.0 * (Gn * (x(3) - N0) - 1.0 / tauP) - kap * eDelay(delayIndex) / x(1) * sin(theta);
                b(i,3) = J - x(3) / tauS - Gn * (x(3) - N0) * x(1) * x(1);
            end
        end
        break
        for i = 1 : M
            a(i) = a(i) + h * (b(1,i) + 2.0 * b(2,i) + 2.0 * b(3,i) + b(4,i)) / 6.0;
        end
        eDelay(delayIndex) = a(1);
        phiDelay(delayIndex) = a(2);
        delayIndex = rem(delayIndex + 1, delayNum);
        if delayIndex == 0
            delayIndex = 1;
        end
    end
    plotCnt = 0;
    prev = 0.0;
    now = 0.0;
    next = 0.0;
    
    for i = 1 : n
        prev = now;
        now = next;
        next = a(1) * a(1);
        if prev <= now && now >= next && i > 2
            o1 = kap * 1e-9;
            o2 = now * 1e-20;
            arrI(length(arrI)+1) = o1;
            arrA(length(arrA)+1) = o2;
            plotCnt = plotCnt + 1;
        	if plotCnt > plotMax 
				break;
            end
        end
        
        t = h * (trans + i);
        x = zeros(1, M);
        b = zeros(4, M);
        theta = rem(phaseShift + a(2) - phiDelay(delayIndex), 2.0 * pi);
        for i = 1 : 4
            for j = 1 : M
                if i == 1
                    x(j) = a(j);
                elseif i == 2
                    x(j) = a(j) + h * b(1,j) / 2.0;
                elseif i == 3
                    x(j) = a(j) + h * b(2,j) / 2.0;
                elseif i == 4
                    x(j) = a(j) + h * b(3,j) / 2.0;
                end
                b(i,1) = 1.0 / 2.0 * (Gn * (x(3) - N0) - 1.0 / tauP) * x(1)+ kap * eDelay(delayIndex) * cos(theta);
                b(i,2) = alpha / 2.0 * (Gn * (x(3) - N0) - 1.0 / tauP) - kap * eDelay(delayIndex) / x(1) * sin(theta);
                b(i,3) = J - x(3) / tauS - Gn * (x(3) - N0) * x(1) * x(1);
            end
        for i = 1 : M
            a(i) = a(i) + h * (b(1,i) + 2.0 * b(2,i) + 2.0 * b(3,i) + b(4,i)) / 6.0;
        end
        eDelay(delayIndex) = a(1);
        phiDelay(delayIndex) = a(2);
        delayIndex = rem(delayIndex + 1, delayNum);
        if delayIndex == 0
            delayIndex = 1;
        end
        end
    end
end
plot(arrI, arrA, '.');