%% Synthetic microstructure generation: Cubic
% We will distribute NxNxN grid points in a cube shaped specimen of size
% LxLxL microns. Then we will create M grains in it using Voronoi
% tessellation. The output will be written to a text file with name
% ms_file_name ready to be sent to the program for simulating virtual
% diffraction patterns. Copy the text file to software-dev account in a
% folder with heXRD detector, material and config files and run
% $> fwd_model config.yml
%
% MTEX and export_fig packages are required. 
%%
L = 1000;                             % Sample dimensions in microns
N = 100;                              % NxNxN uniform grid will be created from the specimen
M = 50;                               % Number of grains
sigma_o = 1.0E-2;                     % Amount of orientation spread (mosaicity) to add to each grain
sigma_e = 5.0E-3;                     % Amount of deformation spread to add to each grain
ms_file_name = 'ms-synth-cubic.csv';  % Name of the file to which the output is written. This file will be the input for forward modeling code.
% Do we want to generate random grains or use previously saved data?
generate_random_grains = 1;           % 0 = Read grain COM, orientation, strain from a  file.
                                      % 1 = Generate random grains
%%
if(generate_random_grains)
    com    = (rand(M, 3) - 0.5) * L;  % Random grain centroid coordinates
    quat   = randq(M);                % Random quaternions for grain orientations
    quat   = squeeze(double(quat));   % Convert from quat MTEX object to normal matrix
    strain = rand(M, 6)/100.0;        % Small random strain at each grain
else
    load NiTi-A_tess;                 % Load a MAT file with com, quat strain
end
% Create a tessellation
tess = zeros(N, N, N);                % A 3D matrix containing grain ID for each point
tess_coords = zeros(N*N*N, 3);        % A 3D matrix containing coordinates of each point
for ii = 1:N
    for jj = 1:N
        for kk = 1:N
            tess_coords(N*N*(ii - 1) + N*(jj - 1) + kk, :) = [ii/N*L - L/2.0 jj/N*L - L/2.0 kk/N*L - L/2.0];
        end
    end
end

% Calculate distance between each material point and grain centroids
dist = pdist2(tess_coords, com);

for ii = 1:N
    for jj = 1:N
        for kk = 1:N
            % Voronoi tessellation is where each material point is assigned
            % to a grain whose COM is closest to that material point
            [~, tess(ii, jj, kk)] = min(dist(N*N*(ii - 1) + N*(jj - 1) + kk, :));
        end
    end
end
%%
% File to write microstructural data
f = fopen(ms_file_name, 'w');

for ii = 1:N
    for jj = 1:N
        for kk = 1:N
            % Write microstructure data in the column format that fwdmodel
            % python code reads.
            % Strain is really right stretch tensor.
            % I am adding a small noise to orientations and strains. Hence
            % the rand() calls.
            fprintf(f, '%12.4f, %12.4f, %12.4f, %s, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f\n', ...
                ii/N*L - L/2.0, jj/N*L - L/2.0, kk/N*L - L/2.0, 'NiTi_cubic', ...
                quat(tess(ii, jj, kk), 1) + rand*sigma_o, quat(tess(ii, jj, kk), 2) + rand*sigma_o, quat(tess(ii, jj, kk), 3) + rand*sigma_o, quat(tess(ii, jj, kk), 4) + rand*sigma_o, ...
                1.0 + strain(tess(ii, jj, kk), 1) + rand*sigma_e, 1.0 + strain(tess(ii, jj, kk), 2) + rand*sigma_e, 1.0 + strain(tess(ii, jj, kk), 3) + rand*sigma_e, ...
                strain(tess(ii, jj, kk), 4) + rand*sigma_e, strain(tess(ii, jj, kk), 5) + rand*sigma_e, strain(tess(ii, jj, kk), 6) + rand*sigma_e);
        end
    end
end
% Don't forget to close files.
fclose(f);
%%
% Plot COM, orientations etc? Useful for debugging.
plot_figures = 1;
%
if(plot_figures)
    % IPF with orientations
    q = quaternion(quat');
    o = orientation(q, crystalSymmetry('cubic'), specimenSymmetry('1'));
    plotIPDF(o, yvector, 'xAxisDirection','east', 'MarkerSize', 8);
    export_fig 'NiTi-A_orientations' -png -r100
    % Scatter plot with tessellation
    figure;
    scatter3(tess_coords(:, 1), tess_coords(:, 2), tess_coords(:, 3), 4, tess(:), 'filled')
    axis equal; axis vis3d;
    export_fig 'NiTi-A_tess' -png -r100
end
