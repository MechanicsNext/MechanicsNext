% Rate independent CP law with hardening
e = 0; % strain
s = 0; % stress
mu = 120e9; % stiffness
gammadot = @(x,y) (0.02*(x/y)^10)*sign(x); % slip rate (Peirce et al)
taucrss0 = 200e6; % Initial tau_CRSS
taucrss = taucrss0;

e_p = 0; % plastic strain
delta_t = 0.0001;
% Loop over 1000 timesteps
for i=1:1000
    e(end+1) = e(end)+delta_t; % Impose 0.0001 strain increment
    s_trial = mu*(e(end) - e_p(end)); % Trial stress using approx elastic deformation
    gammadot_trial = gammadot(s_trial, taucrss); % Calculate slip rate
    e_p(end+1) = e_p(end) + gammadot_trial*delta_t; % Calculate new plastic strain
    taucrss = taucrss + taucrss0*abs(gammadot_trial)*delta_t; % Update hardness
    s(end+1) = mu*(e(end) - e_p(end)); % Calculate new stress using elastic constitutive law
end
% Plot
plot(e, s/1e6, 'Color', 'b', 'LineWidth', 2);
xlabel('strain');
ylabel('Stress (MPa)');