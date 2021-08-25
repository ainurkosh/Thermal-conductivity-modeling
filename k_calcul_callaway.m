% This is code V.1. for numerical calculation of thermal conductivity of 
% SHI irradiated non-metals
clear 

% =========================== parameters ==================================
kB = 1.38064852e-23;            % Boltzman constant in J/K
hr = 6.62607015e-34 / (2 * pi); % reduced Plack constant in J*s
%f  = 10e+12;                    % frequency in Hz
f = 0.8e+12;
w =  2 * pi * f;                % frequency Hz of phonons in LiF crystals 
                                % based on Lindsay 2016
TD = 735;                       % Debye temperature of LiF crystal taken 
                                % from "Atomic vibrations" chapter of 
                                % "Quantum theory of the solid state" by Kantorobich
M = 25.939;                     % Molecular mass LiF in g/mol
g = 1.5;                        % Gruneisen constant of LiF, a measure of anharmonicity
a = 4.03e-10;                   % a^3 atomic volume in m

T = [1 : 100 110 : 10 : 1000];
%T = 300;

% estimated
v = 7000;                       % in the range of 4000 - 7000 m/s based on 
                                % Wright et al. 2005
%dM = 2*M * 1e-1;                % temporary value

L = 2e-3;

% load DFT simulation data
k_num = load('numdata.txt');
% k_num_4A = k_num(:,1);
% k_num_7A = k_num(:,2);
k_num_8A = k_num(:,3);

% load Berman experimental simulation data
k_Ber = load('LiF_k_Berman.txt');
T_Ber = k_Ber(:,1);
k_Ber = k_Ber(:,2)* 1e+2;

% load Pohl experimental simulation data
k_Pohl = load('LiF_k_Pohl.txt');
T_Pohl = k_Pohl(:,1);
k_Pohl = k_Pohl(:,2)* 1e+2;

% _________________________________________________________________________
%                               CALCULATION
% _________________________________________________________________________
% tuning parameters

A = 1e-44;      % pure, unirradiated LiF
%A = 1.1e-44;     % point defects in irradiated LiF starting from 1.1e-44
B = 1.35e-22;
C = 3.83e+13;
dw = 1;
C1 = hr^2 * w^4 * kB^(-1) * T.^(-2);
%C2 = A * w^4 + B * exp(-50./T) .* T.^3 * w^2 + v * L^(-1); % Pohl mentions
%exp(-50./T) term, but it does not improve.
C2 = A * w^4 + B * T.^3 * w^2 + v * L^(-1);
C3 = exp(hr * w ./ (kB .* T));
C4 = (exp(hr * w ./ (kB .* T)) - 1).^2;
k = C / (2 * pi^2 * v) * C1 ./ C2 .* C3 ./ C4 * dw;
display(max(k))
display(find(k == max(k)))
display(k((T == 300)))


figure
lw = 1.3;
%plot( T, k, 'LineWidth', lw)

plot(  T, k, 'b', T, k_num_8A, 'r--',  'LineWidth', lw)
hold on
plot(T_Ber, k_Ber,'s','MarkerSize',5, 'color',[0.0 0.5 0.5],  'LineWidth', lw)
plot(T_Pohl, k_Pohl,'s','MarkerSize',5, 'color',[0.2 0.0 0.7],  'LineWidth', lw)
hold off
legend('k_L_i_F Callaway', 'k_L_i_F DFT', 'experimental Berman 1956', ...
    'experimental Pohl 1960')%, 'experimental Han 1996', 'experimental 9.6% 7Li  Lindsay 2016', ...

xlabel('temperature, K')
ylabel('thermal conductivity')
axis([0 100 0 2850])