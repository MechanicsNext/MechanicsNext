% Macro plasticity model with isotropic hardening
% Harshad Paranjape, March 2016
% hparanja@mines.edu

% Material properties
E               = 120E3;                      % Stiffness (MPa)
H               = 500;                        % Hardening coefficient (MPa)
stress_y0       = 300;                        % Initial yield stress (MPa)
% Simulation parameters
nsteps          = 10000;                      % Total number of time steps
strain_max      = 0.1;                        % Max strain imposed
% Outputs
strain          = zeros(nsteps, 1);           % strain @ each step
strain_plastic  = zeros(nsteps, 1);           % plastic strain @ each step
stress_trial    = zeros(nsteps, 1);           % trial stress @ each step
stress          = zeros(nsteps, 1);           % stress @ each step
stress_y        = zeros(nsteps, 1);           % yield stress @ each step

strain_increment = strain_max/nsteps;
stress_y(1) = stress_y0;

for ii = 1:nsteps-1
    strain(ii + 1) = strain(ii) + strain_increment;
    % Elastic predictor
    stress_trial(ii + 1) = stress(ii) + E*strain_increment;
    % Yield function
    yield_function = stress_trial(ii + 1) - stress_y(ii);
    % Check for yield
    if(yield_function < 0)
        % Elastic
        strain_plastic(ii + 1) = strain_plastic(ii);
        stress_y(ii + 1) = stress_y(ii);
        stress(ii + 1) = stress_trial(ii + 1);
    else
        % Plastic
        delta_lambda = (stress_trial(ii + 1) - stress_y(ii)) / (E + H);                    % Calculate plastic strain increment
        strain_plastic(ii + 1) = strain_plastic(ii) + delta_lambda;                        % Increment plastic strain
        stress(ii + 1) = stress_trial(ii + 1) - E*delta_lambda*sign(stress_trial(ii + 1)); % Plastic corrector
        stress_y(ii + 1) = stress_y(ii) + H*delta_lambda;                                  % Increment yield stress
    end
end

% Plot stress-strain curve
plot(strain, stress, 'b', 'LineWidth', 2);
xlabel('Strain');
ylabel('Stress (MPa)');
set(gca, 'FontSize', 16);
