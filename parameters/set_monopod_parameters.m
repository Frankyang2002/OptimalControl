%% Set monopod parameters
% Init
clear; clc; close all;

%% Elhardt: Table 6.3 & Grimminger Figure 4
% Base parameters
m_B = 0.149;  % Mass [kg]
l_B = 0.160;  % Length [m]
w_B = 0.05;   % Width [m]
I_B = m_B / 12 * (l_B^2 + w_B^2);  % Inertia [kg * m^2]

% Upper leg parameters
m_UL = 0.149;
l_UL = 0.160;
w_UL = 0.04;
I_UL = m_UL / 12 * (l_UL^2 + w_UL^2);

% Lower leg parameters
m_LL = 0.031;
l_LL = 0.160;
w_LL = 0.03;
I_LL = m_LL / 12 * (l_LL^2 + w_LL^2);

%% SLIP-model
% Spring constants
K_l = 300;
K_phi = 1;

% Damper constants
D_l = 6;
D_phi = 1;

% Desired resting spring pose
q_rest = [0; deg2rad(-135); deg2rad(90)];

%% Hopping height control
K_p = 0.1;

%% Enable motor torque limitation
enable_trq_limit = false;
trq_max = 2.7; % see https://arxiv.org/pdf/1910.00093.pdf

%% Save file
if isfile('monopod_parameters.mat')
    save('monopod_parameters.mat')
elseif isfile(fullfile(pwd, 'parameters', 'monopod_parameters.mat'))
    save(fullfile(pwd, 'parameters', 'monopod_parameters.mat'))
end

%% Update parameter changes in dynamics functions
derive_dynamics