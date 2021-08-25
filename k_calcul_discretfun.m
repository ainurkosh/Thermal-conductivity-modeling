% This is code V.1. for estimated calculaton of thermal conductivity of Bi 
% irradiated LiF ceramics. The modeling is based on Klemens theory

% =========================== parameters ==================================

% k - thermal conductivity

% C - contribution of lattice waves of frequency f (f frequency dependent)

% v - grop velocity of the phonons of LiF crystal

% l - mean free pass (f frequency dependent)

% P - reciprocal mean free path or scattering probability per unit path length 

% Pi - scattering probability due to phonon-phonon interactions caused by
% crystal anharmonicity

% Pp - probability of scattering due to point defects

% Px - probability of scattering by extended imperfections such as colloids
% or voids, large compared to the phonon wavelength

% PL - probability of scattering by grain bundaries or extended bpundaries 

f =  10e+12; % frequency (GHz) of phonons in LiF crystals based on Lindsay 2016
t = [1 : 100 110 : 10 : 1000];   % temperature in K
tD = 735; % Debye temperature of LiF crystal taken from "Atomic vibrations" 
          % chapter of "Quantum theory of the solid state" by Kantorobich
          % This value is 620 in "Thermal conductivity and lattice virational modes"
          % by Klemens 1958
g = 1.5; % Gruneisen constant of LiF, a measure of anharmonicity
mu = 55.14e+9; % shear modulus in Pa
a = 4.03e-10; % a^3 atomic volume in m
kB = 1.38064852e-23; % Boltzman constant in J/K
h = 6.62607015e-34; % Plack constant in J*s
M = 25.939; % Molecular mass LiF in g/mol

% estimated parameters
v = 4000; % in the range of 4000 - 7000 m/s based on Wright et al. 2005
C = 4.4e+5; % the contribution of lattice waves (taken to match values of k thermal conductivity) 
c = 1e-3; % the point defect concentration (per atom)
dM = 2*M * 1e-1; % temporary value

% load molecular dynamics simulation data
k_num = load('numdata.txt');
k_num_4A = k_num(:,1);
k_num_7A = k_num(:,2);
k_num_8A = k_num(:,3);

% =========================================================================
% scattering due to phonon-phonon interaction at hight themperatures
t0 = mu * a^3 / kB;
fD = kB * tD / h;
Pi = 4 * pi * g^2 * (t/t0) / v / fD * f^2;


% scattering due to point defect
% A = c * (dM / M)^2 * 4 *pi^3 * a^3 / v^4;
% Pp = zeros(length(A), length(f));
% for i = 1 : length(A)
%     Pp(i, :) = A(i) * f.^4;
% end
Pp = 0;

% scattering by extended defects
% nx = 1e+8; % number of extended defects
% sigma = 1e-6; % cross section (provided that the wavelength is small 
%               % compared to the diameter sqrt(sigma))
% Px = nx * sigma;
Px = 0; % extended defects should be of suffcient concentration in order 
        % to contribute to thermal conductivity degradation, which is not 
        % evidenced in the literature

% scattering by grain and external boundaries
PL = 0;     % strongly decreses thermal conductivity, has especially strong 
            % contribution at lower concentrations of point defects
% PL = 5e+7;
df = 1;
P = Pi;  % + Pp + Px + PL;
l = 1./P; 
k3 = 1/3 * C * v * l * df;

%__________________________________________________________________________
%                       Umklapp processes
%__________________________________________________________________________
% thermal conductivity calculation at umklapp processes
% scattering due to phonon-phonon interaction
% Case A8 (%fitting to A8 LiF mol.dyn. simulation data)
B2 = 2e+28;
alpha = 2.44; % around unit according to Klemens
x = h * f / kB ./ t.^3;
%PU = B2 * 8 * g^2 * kB * tD / M .* x .* exp(- tD ./ ( alpha * t) );
PU = B2 * 8 * g^2 * kB * tD / M * exp(- tD ./ ( alpha * (t + 29)) );

df = 1;
P = PU;  % + Pp + Px + PL;
l = 1./P; 
k2 = 1/3 * C * v * l * df;

%__________________________________________________________________________
% thermal conductivity calculation at lowest temperatures
% scattering due to phonon-phonon interaction
A2 = 2.7e+8;
P0 = A2 ./ t.^3;

df = 1;
P = P0;  % + Pp + Px + PL;
l = 1./P; 
k1 = 1/3 *  C * v * l * df;
% 
figure
plot( t, k1, 'b', t(1:105), k2(1:105), 'm', t(105:end), k3(105:end),'r', t,  k_num_8A, 'k--','LineWidth', 1.5 )
title('k modeling & DFT data fit')
xlabel('temperature (K)')
ylabel('k (W/mK)')
legend('boundary scattering', 'umklapp processes', ...
     'high T U-processes', 'DFT simulation')
axis([0 1000 0 3e+3])
%--------------------------------------------------------------------------
%figure
% plot(t, k_num_4A, 'g--',t, k_num_7A, 'm--', t,  k_num_8A, 'b--', 'LineWidth', 2 )