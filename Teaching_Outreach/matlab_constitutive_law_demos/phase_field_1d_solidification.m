%% A Phase Field Based Model for Solidification
%  Harshad Paranjape, Peter M Anderson, MSE 6756 The Ohio State University
%
% This program models solidification of a melt in a 1D mould.
% The output of the program is the 1D temperature distribution, from 
% position x = x1 to x2, from time t = t1 to t2. The x domain is divided into 
% nx intervals and the t domain into nt intervals.
% Initial and boundary conditions -
% At t = t1, the (initial) temperature distribution is T(x)=T0.
% At x = x1 (in m), the Newtonian boundary condition Jx=h(Tmold-T) holds.
% At x = x2 (in m), the symmetric boundary condition Jx=0 holds.
%
% Outputs
% Temperature, T = T(t, x)
% T is not stored at all timesteps. In stead it is stored for nsteps number
% of equally spaced timesteps.
% Casting parameters are:
% Tmold = Temperature of the mold, assumed to be constant (C)
% h = the Newtonian heat transfer coefficient (W/m^2-K)
% dhf = heat of solidification (J/kg)
% c = heat capacity (J/kg-T)
% k = thermal conductivity (W/m-T)
% r = density (kg/m^3)
% The above properties are assumed to be equal for the liquid and
% solid phases.
%
% Define thermal parameters
x1 = 0;
x2 = 0.025;
nx = 50;
t1 = 0;
t2 = 500;
nt = 500000;
T0 = 1450;
nstore = 50;
TL = 1356;
Tmold = 600;
h = 1000;
dhf = 2e5;
c = 385;
k = 400;
rho = 8.96e3;
% Define phase field parameters
pf_Q = 1e5;
pf_gl = 0;
pf_gs = 2e5;
pf_k = 2.5e-2;
pf_M = 3e-6;
% Calculate the thermal diffusivity (a)
a = k/(rho*c);
hk = h/k;
% We will use explicit time integration for solving the thermal
% conductivity equation. Stability of the explicit solution is dependant on
% the chosen time and spatial step.
% Check the stability condition (Make sure dt is sufficiently small)
dt = (t2 - t1)/nt;
dx = (x2 - x1)/nx;
dtmax = 0.5/a*((x2 - x1)/nx)^2;
if(dt > dtmax)
    disp(['Time step used ' num2str(dt) ' is larger than the ' ...
        'stable time step ' num2str(dtmax) ...
        '. Solution will be unstable. Increase the number of '...
        'time intervals (nt).'])
end
% Implementation of algorithm for solving heat conduction equation (1D)
% Tout stores the output of the code. Columns store the temperature at
% given X but different time instances. Rows store temperature at different
% X but given time.
T = zeros(nstore+1, nx+1);
T_previous = T0*ones(nx+3, 1);
T_current = zeros(nx+3, 1);
T(1, :) = T_previous(2:end-1);
% Order parameter (phi = 0 in solid state, phi = 1 in liquid state)
pf_phi_store = ones(nstore+1, nx+1);
pf_fdrive_store = zeros(nstore+1, nx+1);
pf_dphi_store = zeros(nstore+1, nx+1);
pf_grad_store = zeros(nstore+1, nx+1);
pf_phi = pf_phi_store(1, :)';
% xsteps stores the x-locations at which temperature is calculated.
xsteps = x1:(x2-x1)/nx:x2;
% tsteps stores time steps at which temperature is calculated.
tsteps = zeros(nstore+1, 1);
%
ntsteps = floor(nt/nstore);
% Heat capacity
c2 = a*dt/dx^2;
%
counter = 2;
for i=1:nt
    % Update the order parameter from PF kinetic law    
    % Define free energy function
    pf_h = (pf_phi.^2).*(3.0 - 2.0*pf_phi);
    pf_f = (pf_phi.^2).*(1 - pf_phi).^2;
    % Define driving force as the derivative of free energy wrt the order
    % parameter
    pf_h_prime = 6*pf_phi - 6*pf_phi.^2;
    pf_f_prime = 2*pf_phi - 6*pf_phi.^2 + 4*pf_phi.^3;    
    % Bulk free energy of the solid-liquid system
    pf_g = ((TL - T_previous(2:end-1))/TL).*(pf_h*pf_gs + (1-pf_h)*pf_gl) + pf_Q*pf_f;    
    % Define the gradient term
    pf_grad_term = pf_k*del2(pf_phi, (x2 - x1)/nx);
    % Since the first node does not have a neighbor, the gradient term is
    % incorrectly calculated there. Set it to zero.
    pf_grad_term(1) = 0;
    % Get the total driving force for solidification = bulk contribution +
    % gradient term
    pf_fdrive = ((TL - T_previous(2:end-1))/TL).*(pf_h_prime*(pf_gl - pf_gs)) + pf_Q*pf_f_prime + pf_grad_term;
    % d(phi)/dt = -M*fdrive
    % Obtain the increment in the order parameter
    pf_dphi = pf_M*pf_fdrive;
    % Add noise at the first node. This simulates nucleation of the solid
    % phase
    %pf_dphi(1) = pf_dphi(1) + (rand(1, 1) - 0.5)/1.0e3;
    pf_dphi = pf_dphi + (rand(size(pf_dphi)) - 0.5)/1.0e6;
    % Get the new order parameter at the current time step
    pf_phi = pf_phi + dt*pf_dphi;
    % Check for bounds on the order paramer. 0 <= phi <= 1
    pf_phi(pf_phi(:) > 1) = 1;
    pf_phi(pf_phi(:) < 0) = 0;
    % Update the temperature at each node to account for the energy
    % released or absorbed by the order parameter change
    for j=2:nx+2
        T_previous(j) = T_previous(j) - dt*pf_dphi(j-1)*dhf/c;
    end
    % Enforce the zero flux BC at the right end by setting the temperature 
    % at the two right most nodes equal.
    T_previous(end) = T_previous(end-1);
    % Calculate temperature of left most dummy node to take into account
    % convective cooling at the mold surface. This way we can uniformly use
    % conduction equation at all the real nodes to calculate temperature
    % there.
    T_previous(1) = T_previous(2) + hk*dx*(Tmold - T_previous(2));    
    for j=2:nx+2
        % Find the temperature at the nodes by solving the conduction heat
        % equation using the finite difference method.
        T_current(j) = T_previous(j) + c2*(T_previous(j+1) - ...
            2*T_previous(j) + T_previous(j-1));
    end
    % Store temperature, order parameter, driving force etc at for nstore
    % number of iterations
    if(mod(i, ntsteps) == 0)
        T(counter, :) = T_current(2:end-1);
        pf_phi_store(counter, :) = pf_phi';
        pf_dphi_store(counter, :) = dt*pf_dphi';
        pf_fdrive_store(counter, :) = pf_fdrive';
        pf_grad_store(counter, :) = pf_grad_term';
        tsteps(counter) = t1 + (counter - 1)*(t2 - t1)/nstore;
        % Display progress
        disp(['Saved output number ' num2str(counter) ' of ' num2str(nstore+1)])
        counter = counter + 1;        
    end
    T_previous = T_current;
end
% Visualizations
% Plot temperature at the last time instant across the length of the
% sample.
figure;
plot(xsteps, T(end, :));
xlabel('Position (m)');
ylabel('Temperature (K)');
title('Temperature variation with position @ last time step');
% Plot temperature variation with time at the right end of the sample
figure;
plot(tsteps, T(:, end));
xlabel('Time (s)');
ylabel('Temperature (K)');
title('Temperature variation with time @ right end of the sample');
% Plot order param variation with time at the right end of the sample
figure;
plot(tsteps, pf_phi_store(:, end));
xlabel('Time (s)');
ylabel('Order parameter');
title('Order parameter variation with time @ right end of the sample');
figure;
pcolor(pf_phi_store);
title('Order parameter');
% Additional plots for debugging
% figure;
% pcolor(pf_dphi_store);
% title('Order parameter increment');
% figure;
% pcolor(pf_fdrive_store);
% title('Driving force');
% figure;
% pcolor(pf_grad_store);
% title('Gradient term');