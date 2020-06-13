clc;clear;
%% Parameters calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%physical constants%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kb = 8.617e-5; % eV/K
q  = 1.602e-19; % C
m0=9.11e-31; %kg
meff = 2.07e-32; % kg
h = 4.13e-15; % eV.s
hba=6.582e-16; %eV.s
c=2.9979e10; %cm/s
epsilon0=8.854e-14; %C^2/(J.cm)
%%%%%%%%%%%%%%%%%%%%%%%%%%parameters for material%%%%%%%%%%%%%%%%%%%%%%%%%%
T  = 295;  % K % temperature %290,300
E_WL = 0.97; % eV % energy level of WL %0.97
E_ES = 0.87; % eV % energy level of ES %0.87
E_GS = 0.82; % eV % energy level of GS %
wa_WL=E_WL/hba;  %% angular frequency
wa_ES=E_ES/hba;
wa_GS=E_GS/hba;
nr = 3.5;
R1 = 0.32;  % 0.32
R2 = R1;
alpha_in = 6; % cm^-1 % internal loss %6,10
Beta = 1e-4; % spontaneous emission coefficient %1e-3
Gamma = 0.06; % optical confinement factor %
tau_WL_spon = 500e-12;  % s % spontaneous emission time for WL %
tau_ES_spon = 500e-12;  % s % spontaneous emission time for ES %
tau_GS_spon = 1200e-12; % s % spontaneous emission time for GS %
T_D= 1e-13; %% dephasing time unit=s 0.1ps
h_WL=5e-7;  % cm % WL thickness%   origianl =1e-7
h_QD=5e-7;  % cm % qd thickness per layer%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 5e-2; % cm % cavity length %
W = 4e-4;  % cm % cavity width %
n_layer = 5; % number of QD layers %
ND = 10e10;  % cm^-2 % QD surface density %
eta_tau=0.25; % normally 0.25
tau_WL_ES = eta_tau*25.1e-12; %s  12.6 ps, 6.3 ps
tau_ES_GS = eta_tau*11.6e-12; %s 5.8 ps, 2.9 ps
a_GS =5e-15; %cm^2% gain coefficent for GS %  5e-15
a_ES =2.*a_GS; %cm^2% gain coefficient for ES % fitting parameter 2
a_WL=0.5.*a_GS; % cm2  differential gain for WL 0.5
epsilon_V= 2e-16; %2e-16;% cm^3  gain compression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%Parameters calculations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% For GS lasing
F_WL=wa_GS./wa_WL.*(wa_WL-wa_GS).*T_D/(1+(wa_WL-wa_GS)^2*T_D^2); % coefficient of refractive index conrtibution to GS
F_ES=wa_GS./wa_ES.*(wa_ES-wa_GS).*T_D/(1+(wa_ES-wa_GS)^2*T_D^2); % coefficient of refractive index conrtibution to GS
tau_p =(c.*(alpha_in + log(1./(R1.*R2))./(2.*L))./nr).^-1; % photon lifetime9.9940e-012  ,5.8ps
S = L*W; % cm^2 % Surface area %3.3e-5
Nb = ND*S*n_layer;% QD number % 8.25e6
V_QD=S*h_QD*n_layer; % volume for active layer%8.25e-11
V_WL=S*h_WL; % volume for Wetting layer%3.3e-12
rho_WL1=meff.*Kb.*T.*S./(pi.*(hba.^2));%% unit is Kg*(cm/s)^2/eV
rho_WL = rho_WL1.*(1e-4)./q; %unit is Kg*(m/s)^2/J, or unit is 1  for denegeracy of WL
epsilon_GS = epsilon_V.*Gamma./V_QD; % gain compression factor for gS%
epsilon_ES = 0.*epsilon_GS; % gain compression factor for eS%1e-7
tau_GS_ES = 1/2.*tau_ES_GS.*exp((E_ES - E_GS)./(Kb.*T)); % Carrier escape time from GS to ES
tau_ES_WL = 4.*tau_WL_ES.*exp((E_WL - E_ES)./(Kb.*T)).*Nb./rho_WL;  % Carrier escape time from ES to WL;
%%%%%%%%%%%%%%%%% Feedback parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_LEF= 1;
tau_nr=10e-9; %s,non-radiative recombination lifetime Qdot lasers are 1e-8 s to 2e-8 s
tau_in = 2*nr*L/c;% internal round trip delay
R_ref = 1e-8; %1e-6,1e-4,1e-2%% feedback fraction 
k_c = (1-R1).*sqrt(R_ref)./sqrt(R1); % feedback strength
L_ex = 0.35; % cm; 0.1, 100;
n_ex = 1.5;% refractive index in the fiber;
tau_ex = 2*L_ex*n_ex/c;% roundtrip time in the external cavity;
phaseShift= 2*pi.*tau_ex.*E_GS./h; % feedback phase, 0.0098...0.98*2*pi; then the phase is always about 2*pi; if Lex=0.15 cm, 50 cm, then the phase is about pi;
phaseShift = rem(phaseShift,2*pi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Runge Kutta 4 slover
M = 6;
h = 5.0e-12;
transient = 5000.0e-9; 
tMax = 20.0e-9; 
div = 10; 
trans = floor(transient / h);
n = floor(tMax / h);

%a - initial condition - 待定 等算的时候我给你说
% nr - 10e-9, 2*Ith(52mA) = 104mA
a = [7.231082211578302e+05, 1.215813271872004e+07, 1.496669230548948e+07, 0, 1.465414777198156e+05, 0];
%a = [7.231082211578302e+05, 1.215813271872004e+07, 1.496669230548948e+07, 0, 1.465414777198156e+05, 0];


% initial phase shift
delayNum = floor(tau_ex / h);
delayIndex = 1;
DELAY_MAX = 100000;
eDelay = zeros(1, DELAY_MAX);
phiDelay = zeros(1, DELAY_MAX);

for i = 1 : DELAY_MAX
    eDelay(i) = a(5);
    phiDelay(i) = a(6);
end

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
            %x(1) - Nwl
            %x(2) - Nes
            %x(3) - Ngs
            %x(4) - Ses = 0
            %x(5) - Sgs
            %x(5) - phi
            I=0.104;
            f_ES = 1 - x(2)/(4*Nb);  % non-occupation probabilityfor ES
            f_GS = 1 - x(3)/(2*Nb);  % non-occupation probabilityfor GS
            
            b(i,1) = I/q + x(2)/tau_ES_WL - x(1)*f_ES/tau_WL_ES - x(1)/tau_WL_spon - x(1)/tau_nr;
            b(i,2) = x(1)*f_ES/tau_WL_ES + x(3)*f_ES/tau_GS_ES - x(2)/tau_ES_WL - x(2)*f_GS/tau_ES_GS - x(2)/tau_ES_spon - x(2)/tau_nr;
            b(i,3) = x(2)*f_GS/tau_ES_GS - x(3)*f_ES./tau_GS_ES - x(3)/tau_GS_spon - x(3)/tau_nr - Nb/V_QD*Gamma*c/(nr)*a_GS*(2*x(3)/(2*Nb) - 1)*x(5)/(1 + epsilon_GS*x(5));
            b(i,4) = 0;
            b(i,5) = Nb*Gamma*c/(nr)*a_GS/V_QD*(2*x(3)/(2*Nb) - 1)*x(5)/(1 + epsilon_GS*x(5)) - x(5)/tau_p + Beta*x(3)/tau_GS_spon + 2*k_c/tau_in*sqrt(x(5)*eDelay(delayIndex))*cos(theta);  
            b(i,6) = 0.5*alpha_LEF*(Gamma*c/(nr)*a_GS/V_QD*Nb*(2*x(3)/(2*Nb) - 1)/(1 + epsilon_GS*x(5))-1/tau_p) - k_c/tau_in*sqrt(eDelay(delayIndex)/x(5))*sin(theta);
        end
    end
    break
    for i = 1 : M
        a(i) = a(i) + h * (b(1,i) + 2.0 * b(2,i) + 2.0 * b(3,i) + b(4,i)) / 6.0;
    end
    eDelay(delayIndex) = a(5);
    phiDelay(delayIndex) = a(6);
    delayIndex = delayIndex + 1;
end

arrI = [];
arrA = [];

for i = 1 : n
    t = h * (trans + i);
    % output i % div == 0
    o1 = h * i * 1e9;
    o2 = a(5);
    arrI(i) = o1;
    arrA(i) = o2;
    
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
            I=0.104;
            f_ES = 1 - x(2)/(4*Nb);  % non-occupation probabilityfor ES
            f_GS = 1 - x(3)/(2*Nb);  % non-occupation probabilityfor GS
            
            b(i,1) = I/q + x(2)/tau_ES_WL - x(1)*f_ES/tau_WL_ES - x(1)/tau_WL_spon - x(1)/tau_nr;
            b(i,2) = x(1)*f_ES/tau_WL_ES + x(3)*f_ES/tau_GS_ES - x(2)/tau_ES_WL - x(2)*f_GS/tau_ES_GS - x(2)/tau_ES_spon - x(2)/tau_nr;
            b(i,3) = x(2)*f_GS/tau_ES_GS - x(3)*f_ES./tau_GS_ES - x(3)/tau_GS_spon - x(3)/tau_nr - Nb/V_QD*Gamma*c/(nr)*a_GS*(2*x(3)/(2*Nb) - 1)*x(5)/(1 + epsilon_GS*x(5));
            b(i,4) = 0;
            b(i,5) = Nb*Gamma*c/(nr)*a_GS/V_QD*(2*x(3)/(2*Nb) - 1)*x(5)/(1 + epsilon_GS*x(5)) - x(5)/tau_p + Beta*x(3)/tau_GS_spon + 2*k_c/tau_in*sqrt(x(5)*eDelay(delayIndex))*cos(theta);  
            b(i,6) = 0.5*alpha_LEF*(Gamma*c/(nr)*a_GS/V_QD*Nb*(2*x(3)/(2*Nb) - 1)/(1 + epsilon_GS*x(5))-1/tau_p) - k_c/tau_in*sqrt(eDelay(delayIndex)/x(5))*sin(theta);
        end
    for i = 1 : M
        a(i) = a(i) + h * (b(1,i) + 2.0 * b(2,i) + 2.0 * b(3,i) + b(4,i)) / 6.0;
    end
    eDelay(delayIndex) = a(5);
    phiDelay(delayIndex) = a(6);
    delayIndex = delayIndex + 1;
    end
end
%% PLOT AREA
plot(arrI, arrA);