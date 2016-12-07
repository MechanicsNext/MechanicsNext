%% Cellular Automaton Model of Slip Patches--PM Anderson, Harshad Paranjape
%   --MSE756 -- Ohio State University
%   A cube of dimension h is deformed with a uniform shear rate. A cross 
%   sectional area of h x h is discretized into nx by ny patches, each of 
%   which has a critical shear strength to initiate a plastic slip of 
%   magnitude b in that patch. Such events redistribute stress, causing 
%   slip in neighboring regions. The output is analyzed in terms of the 
%   number and location of slip events and the stress redistribution.
%
% Define the number of elements along the x and y directions
nx = 40;
ny = 40;
% Define a matrix of nx x ny elements containing the shear strengths (in MPa)
% containing shear strengths that follow a normal distribution.
savg = 100; % Average
sdev = 5;  % Standard deviation
smax = 200;
Tc = savg + sdev*randn(ny, nx);
% Assign values to bh = b/h;
bh = 2.0E-5;
% The elastic shear modulus mu (in MPa)
mu = 70.0E3;
% Impose a macro shear strain increment on the assembly that is 0.3 of the 
% maximum plastic strain increment, attainable if all the patches slip in 
% one increment.
deltagamma = 0.3*bh;
% If a patch slip, the increment in plastic strain is deltagamma_local =
% b/hz; where hz = h/sqrt(nx*ny). The corresponding shear stress drop is
% deltasigma_local = 2/3*mu*deltagamma_local. Based on the Eshelby solution
% for an ellipsoid that shears by deltagamma_local.
deltagamma_local = bh*sqrt(nx*ny);
% The 2/3 factor varies with Poisson's ratio
deltasigma_local = 2.0/3.0*mu*deltagamma_local;
% Write a subroutine that computes the matrix T containing the shear stress 
% in each patch and a matrix Gamma containing the number of slip events in 
% each patch, after niter iterations. The input is the starting T and Gamma
% matrices and the output is Load = [T:Gamma] after niter. Initial macro
% stress = T(0,0), initial total number of slip events = Gamma(0,0),
% initial increment in total slip events = Gamma(0,1)
% The initial T matrix is a uniform stress equal to the shear strength of 
% the weakest patch minus 2MPa
Tmin = min(Tc(:)) - 2.0;
T = Tmin*ones(ny, nx);
T_0_0 = Tmin;
% The intial number of slip events is zero in all patches
Gamma = zeros(ny, nx);
Gamma_0_0 = 0;
Gamma_0_1 = 0;
% Number of iteration
niter = 120;
%
tau = zeros(niter, 1);    % tau = macro stress
np = zeros(niter, 1);     % np = cumulative number of slip events
npinc = zeros(niter, 1);  % npinc = number of slip events in a single iter
tau(1) = T_0_0;
np(1) = Gamma_0_0;
npinc(1) = Gamma_0_1;
% T_store and Gamma_store are defined to store local stress and slip events
% respectively for all iterations.
% Later you can use these variables (e.g. T_store(:, :, 20)) to plot local
% distributions of those quantities at a given iteration.
T_store = zeros(ny, nx, niter);
Gamma_store = zeros(ny, nx, niter);
h1 = figure;
for i=2:niter
    Tnew = zeros(ny, nx);
    % Calculate increment in shear stress
    if i == 2
        tauinc = (deltagamma - npinc(1)/(nx*ny)*bh)*mu;
    else
        tauinc = (deltagamma - (np(i-1) - np(i-2))/(nx*ny)*bh)*mu;
    end
    tau(i) = tau(i-1) + tauinc;
    np(i) = np(i-1);
    for ii=1:ny
        for jj=1:nx
            T(ii, jj) = T(ii,jj) + tauinc; % The macro shear increment is intially applied to all patches
            Tnew(ii,jj) = 0;
        end
    end
    %
    for ii=1:ny
        for jj=1:nx
            % If the local crit shear is exceeded, then...
            if(T(ii,jj) > Tc(ii,jj))
                np(i) = np(i) + 1;
                Gamma(ii,jj) = Gamma(ii,jj) + 1;
                % Reduce the stress in the slipped patch
                Tnew(ii,jj) = Tnew(ii,jj) - deltasigma_local;
                % Strain hardening - increase the Tc in the grain when it
                % slips
                Tc(ii, jj) = Tc(ii, jj) + 50*deltagamma_local*abs((smax - Tc(ii, jj)));
                % Get coordinates of four neighbors of the patch
                if(ii == 1)
                    im = ny;
                else
                    im = ii - 1;
                end
                if(ii == ny)
                    ip = 1;
                else
                    ip = ii + 1;
                end
                if(jj == 1)
                    jm =nx;
                else
                    jm = jj - 1;
                end
                if(jj == nx)
                    jp = 1;
                else
                    jp = jj + 1;
                end
                % Redistribute stress to neighbors
                Tnew(im,jj) = Tnew(im,jj) + 0.25*deltasigma_local;
                Tnew(ip,jj) = Tnew(ip,jj) + 0.25*deltasigma_local;
                Tnew(ii,jm) = Tnew(ii,jm) + 0.25*deltasigma_local;
                Tnew(ii,jp) = Tnew(ii,jp) + 0.25*deltasigma_local;
            end
        end
    end
    npinc(i) = np(i) - np(i-1);
    T = T + Tnew;
    T_store(:, :, i) = T;
    Gamma_store(:, :, i) = Gamma;
    imagesc(T);
    colorbar; axis square; axis off; caxis([0 140]);
    drawnow;
    T_movie(i) = getframe;
end
T_0_0 = tau(niter);
Gamma_0_0 = np(niter);
Gamma_0_1 = npinc(niter);
%
macrostress = tau;               % Macro stress
macrostrain = np*bh/(nx*ny);     % Macro strain
strain_inc = npinc;              % Strain increment
% Plot of macrostress vs strain
figure;
plot(macrostrain, macrostress);
xlabel('Macro strain');
ylabel('Macro stress');
% Plot of number of slip events in an increment vs increment number
figure;
plot(1:niter, strain_inc)
xlabel('Iterations');
ylabel('Number of slip events');
% Visualize T_store, Gamma_store and Tc for all elements
% figure; imagesc(Tc);
% figure; imagesc(T_store(:, :, 55));
% figure; imagesc(Gamma_store(:, :, 55));
