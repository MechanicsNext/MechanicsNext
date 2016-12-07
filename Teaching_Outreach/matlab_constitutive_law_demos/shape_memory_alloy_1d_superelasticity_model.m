% Simple, Phenomenological, 1D model for superelasticity
%  Based on A. Sadjadpour, K. Bhattacharya, A micromechanics inspired 
%  constitutive model for shape-memory alloys: the one-dimensional case, 
%  Smart Mater. Struct. 16 (2007) S51. doi:10.1088/0964-1726/16/1/S06.
%
% MATLAB demonstration by Harshad Paranjape
% 3/29/2016
%
clear all;
close all;
% Material parameters
E              = 110e3;                         % Young's modulus (MPa)
theta_as       = 253;                           % Transformation temperatures (K)
theta_ms       = 173;
theta_t        = (theta_ms + theta_as)/2.0;     % A-M equilibrium temperature
theta          = 323;                           % Test temperature (K) must be > theta_as
L              = 80;                            % Latent heat (MJ/m^3)
lambda_dot_0   = 20;                            % Kinetic coefficient in the martensite volume evolution equation
d_lambda_crit  = L*(theta_as - theta_ms)/(theta_as + theta_ms); % Critical driving force for phase transformation
% Simulation variables
num_increments = 600000;                        % Number of time steps in the simulation (should be large)
stress         = zeros(num_increments, 1);      % Stress (MPa)
strain         = zeros(num_increments, 1);      % Total strain
strain_max     = 0.09;                          % Maximum imposed strain in the simulation
delta_strain   = 2.0*strain_max/num_increments; % Strain increment
t              = 100.0;                         % total time for the simulation (generally strain_max/t < 10^-3)
delta_t        = 2.0*t/num_increments;
strain_m       = zeros(num_increments, 1);      % Transformation strain
strain_m_t     = 0.08;                          % Transformation strain bound for tension
strain_m_c     = -0.08;                         % Transformation strain bound for compression
d_lambda       = zeros(num_increments, 1);      % Driving force for phase transformation
lambda         = zeros(num_increments, 1);      % Martensite volume fraction
% Loading or unloading flags
is_loading = 1;                                 % 1 = loading in progress
is_unloading = 0;
% Loop over all time steps
strain_m(1) = strain_m_t;
for ii = 2:num_increments
    if(ii > num_increments/2.0 && is_loading == 1 && is_unloading == 0)
        % When we reach halfway in terms of timesteps, swich from loading
        % to unloading by reversing the sign of strain increment
        is_loading = 0;
        is_unloading = 1;
        delta_strain = -1.0*delta_strain;
    end
    % Set strain_m to appropriate value based on the stress
    if(stress(ii - 1) < 0)
        strain_m(ii) = strain_m_c;
    else
        strain_m(ii) = strain_m_t;
    end
    % Elastic predictor
    strain(ii) = strain(ii - 1) + delta_strain;
    stress(ii) = E*(strain(ii) - lambda(ii - 1)*strain_m(ii - 1));
    % Calculate the driving force for phase transformation
    d_lambda(ii) = stress(ii)*strain_m(ii) - L/theta_t*(theta - theta_t);
    % Martensite increment, based on the driving force
    lambda_dot = 0;
    if(d_lambda(ii) > d_lambda_crit && is_loading == 1)
        lambda_dot = lambda_dot_0*(d_lambda(ii) - d_lambda_crit);
    elseif(d_lambda(ii) < -1.0*d_lambda_crit && is_unloading == 1)
        lambda_dot = lambda_dot_0*(d_lambda(ii) + d_lambda_crit);
    end
    %
    lambda(ii) = lambda(ii - 1) + lambda_dot*delta_t;
    % Performs checks on the bounds of lambda
    % Martensite fraction can never be outside [0, 1]
    if(lambda(ii) > 1)
        lambda(ii) = 1;
    elseif(lambda(ii) < 0)
        lambda(ii) = 0;
    end
    % Transformation corrector
    stress(ii) = E*(strain(ii) - lambda(ii)*strain_m(ii));
end
% Plot
figure;
plot(strain, stress, 'b', 'LineWidth', 2);
xlim([0 strain_max]);
xlabel('Strain');
ylabel('Stress (MPa)');
set(gca, 'FontSize', 16);
%
figure;
plot(strain, lambda, 'm', 'LineWidth', 2);
xlim([0 strain_max]);
xlabel('Strain');
ylabel('Martensite fraction');
set(gca, 'FontSize', 16);
