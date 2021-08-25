% This is code V.1. for numerical calculation of thermal conductivity of 
% SHI irradiated non-metals
clear 

% =========================== parameters ==================================
kB = 1.38064852e-23;            % Boltzman constant in J/K
hr = 6.62607015e-34 / (2 * pi); % reduced Plack constant in J*s
f  = 10e+12;                    % frequency in Hz
w =  2 * pi * f;                % frequency Hz of phonons in LiF crystals 
                                % based on Lindsay 2016
TD = 735;                       % Debye temperature of LiF crystal taken 
                                % from "Atomic vibrations" chapter of 
                                % "Quantum theory of the solid state" by Kantorobich
M = 25.939;                     % Molecular mass LiF in g/mol
g = 1.5;                        % Gruneisen constant of LiF, a measure of anharmonicity
a = 4.03e-10;                   % a^3 atomic volume in m
d = 1;                          % dimensions of the sample  

T = [1 : 100 110 : 10 : 1000];

% estimated
v = 7000;                       % in the range of 4000 - 7000 m/s based on 
                                % Wright et al. 2005
dM = 2*M * 1e-1;                % temporary value


% load DFT simulation data
k_num = load('numdata.txt');
k_num_4A = k_num(:,1);
k_num_7A = k_num(:,2);
k_num_8A = k_num(:,3);


% load Berman experimental simulation data
k_Ber = load('LiF_k_Berman.txt');
T_Ber = k_Ber(:,1);
k_Ber = k_Ber(:,2)* 1e+2;

% load Pohl experimental simulation data
k_Pohl = load('LiF_k_Pohl.txt');
T_Pohl = k_Pohl(:,1);
k_Pohl = k_Pohl(:,2)* 1e+2;

% =========================== calculation =================================
% kL - thermal conductivity of the lattice
c1 = kB / (2 * pi^2 * v);       % constant coefficients 
c2 = (kB / hr)^3;

x = hr * w ./ (kB * T);         % dimensionless parameter

% ------------------phonon scattering relaxation time ----------------- 

% =========================================================================
% FITTING PARAMETERS
c = 1e-11;                       % the point defect concentration (per atom)

% coefficients
% v.1.  (Tritt 2004)
a0 = 1e+39;                      
alpha = 1.6;
% % v.2. (Klemens)
% a0 = 1e+27;
% alpha = 1.6;

% scaling
% a2 = 1.93335e+11;     %fit by peak to DFT
a2 = 0.57e+11;
%==========================================================================

% Umklapp scattering
% v.1. according to Tritt
tU = (a0 * hr * g^2 * w^2 * T .* exp(-TD ./ (alpha*T)) / (M * v^2 * TD)).^(-1);
% ***
% % v.2. according to Klemens
% tU = (a0 * 8 * g^2 * kB * TD * v / M * T .* exp(- TD ./ ( alpha * T ))).^(-1);

% scattering due to point defect
A = c * (dM / M)^2 * 4 * pi^3 * a^3 / v^3;
tPD = (A * f^4).^(-1);     % t^(-1) = v * P (probability by Klemens)

% boundary scattering
tB = (v / d)^(-1);

% summation of photon relaxation times
tq = (tU.^(-1) + tPD.^(-1) + tB.^(-1)).^(-1);

% phonon-phonon normal scattering the relaxation rate
a1 = 1e-11;  % constant
tN = 1 ./ (a1 * w * T.^3);   % for LiF and diamond

% combined relaxation time
tc = (tq.^(-1) + tN.^(-1)).^(-1);    

dx = 1;
I1 = tc .* x.^4 .* exp(x) ./ (exp(x) - 1).^2 * dx;
I2 = tc ./ tN .* x.^4 .* exp(x) ./ (exp(x) - 1).^2 * dx;
I3 = tc ./ (tN .* tq) .* x.^4 .* exp(x) ./ (exp(x) - 1).^2 * dx;


k1 = c1 * c2 * T.^3   .*   I1;
k2 = c1 * c2 * T.^3   .*   I2.^2 ./ I3;
k = a2 * (k1 + k2);

% I = tq .* x.^4 .* exp(x) ./ (exp(x) - 1).^2 * dx;
% kL = c1 * c2 * T.^3 .* I;

% display the temperature @ peak
%display(find(kL == max(kL)))
% display (max(k_num_8A))
% display (max(kL))

%figure
lw = 1.3;
plot(  T, k_num_8A, 'r--', T, k, 'b',  'LineWidth', lw)
hold on
plot(T_Ber, k_Ber,'-s','MarkerSize',5, 'color',[0.0 0.5 0.5],  'LineWidth', lw)
plot(T_Pohl, k_Pohl,'s','MarkerSize',5, 'color',[0.2 0.0 0.7],  'LineWidth', lw)

title('LiF T dependent thermal conductivity')
xlabel('temperature, K')
ylabel('thermal conductivity, W/mK')
legend('k_L_i_F DFT', 'k_L_i_F Debye-Klemens', 'experimental Berman 1956', ...
    'experimental Pohl 1960') 
hold off

%plot(T, k, 'k',  'LineWidth', lw)
%title('LiF T dependent thermal conductivity Debye-Klemens')
%xlabel('temperature, K')
%ylabel('thermal conductivity, W/mK')
%axis([0 300 0 3000])