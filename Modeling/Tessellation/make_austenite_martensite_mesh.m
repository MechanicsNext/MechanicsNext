% MAKE AUSTENITE-MARTENSITE MESH
% Harshad Paranjape | mechanicsNext
% hparanja@mines.edu
%
% Generate a 3D mesh consisting of an Austenite region and an adjacent
% region with twinned martensite consisting of two variants (CVs).
%
% INPUTS
% 1. Twin, invariant plane normals, CV volume ratio.
% 2. Two additional planes to superimpose (optional).
% 3. Colors for various planes
% 4. A mesh input (optional).
%
% OUTPUTS
% 1. A graphical representation of the microstructure (with figure handle).
% 2. An array of element centroids (Nx3).
% 3. An array of integers representating the phase at each element (Nx1).
% 4. A microstructure file that can be used in the far-field diffraction
% simulation framework to create a synthetic diffraction pattern.
%
% NOTES:
% 1. At this time meshses with mixed element types 
% (e.g., brick + tetrahedron) are not processed.
%
% PREREQUISITES
% Geom3D: https://www.mathworks.com/matlabcentral/fileexchange/24484-geom3d
%
%% INPUTS
habit_plane_normal_ip = [-0.5624 -0.7644 -0.3152];                         % Normal to the habit plane (m) in global coordinates
twin_plane_normal_ip  = [0.9736 -0.2281 0];                                % Twin plane normal (n) in the global coordinates
CV_volume_ratio_ip    = 0.71;                                              % CV1:CV2 volume ratio
slip_plane_normal     = [];                                                % Slip plane normal in global coordinates (Optional, leave empty to not draw)
slip_direction        = [];                                                % Slip direction in global coordinates (Optional, leave empty to not draw)
foil_plane            = [];                                                % TEM foil plane normal (Optional, leave empty to not draw)
%
austenite_color_ip    = 'g';                                               % Austenite color in visualization
CV1_color_ip          = 'r';                                               % CV1 color in visualization
CV2_color_ip          = 'b';                                               % CV2 color in visualization
slip_plane_color_ip   = [0.6 0.6 0.6];                                     % Slip plane color
slip_plane_opacity    = 0.4;                                               % Slip plane opacity
%
bbox_half_width       = 2;                                                 % Half width of the bounding box. The box will span -width to +width in X, Y, Z
number_of_twins       = 5;                                                 % Number of twin pairs in the martensite region to draw.
%
mesh_data_dir         = 'example_make_austenite_martensite_mesh';          % Directory where node and connectivity data live
node_file_name        = 'nodes.inp';                                       % Node file with 4 columns - node number, X, Y, Z coordinate
connectivity_file_name= 'connectivity.inp';                                % Connectivity file with M + 1 columns - Element number, number of M elements that define the connectivity
elems_input           = [];                                                % Variable defining coordinates of tessellation mesh grid (Mx3 array) (Optional. Leave empty if using files)
%
q_austenite           = [1.0 0 0 0];                                       % Orientation of austenite in quaternion
q_cv_1                = [1.0 0 0 0];                                       % Orientation of CV 1 in quaternion
q_cv_2                = [1.0 0 0 0];                                       % Orientation of CV 2 in quaternion
sigma_o = 1e-5;                                                            % Orientation spread in the tessellation data
sigma_e = 1e-4;                                                            % Strain spread in the tessellation data
%
%% Calculate a graphical representation of the microstructure and assign
% appropriate phases to each element in the mesh.
[h, elems, elems_phase] = generate_habit_twin_slip_planes_3d(austenite_color_ip, CV1_color_ip, CV2_color_ip, ...
                                                      habit_plane_normal_ip, twin_plane_normal_ip, ...
                                                      CV_volume_ratio_ip, slip_plane_normal, slip_direction, ...
                                                      slip_plane_color_ip, slip_plane_opacity, foil_plane, ...
                                                      bbox_half_width, number_of_twins, ...
                                                      mesh_data_dir, node_file_name, connectivity_file_name, ...
                                                      elems_input);
%
%% Print microstructure file
f = fopen('ms-am-cv10-11.csv', 'w');
for ii = 1:size(elems, 1)
    if(elems_phase(ii) == 2)
    fprintf(f, '%14.8f, %14.8f, %14.8f, %s, %14.8f, %14.8f, %14.8f, %14.8f, %14.8f, %14.8f, %14.8f, %14.8f, %14.8f, %14.8f\n', ...
                500 * elems(ii, 1), 500 * elems(ii, 1), 500 * elems(ii, 1), 'NiTi_cubic', ...
                q_austenite(1) + rand*sigma_o, q_austenite(2) + rand*sigma_o, q_austenite(3) + rand*sigma_o, q_austenite(4) + rand*sigma_o, ...
                1.0 + rand*sigma_e, 1.0 + rand*sigma_e, 1.0 + rand*sigma_e, ...
                rand*sigma_e, rand*sigma_e, rand*sigma_e);
    elseif(elems_phase(ii) == 1)
        fprintf(f, '%14.8f, %14.8f, %14.8f, %s, %14.8f, %14.8f, %14.8f, %14.8f, %14.8f, %14.8f, %14.8f, %14.8f, %14.8f, %14.8f\n', ...
                500 * elems(ii, 1), 500 * elems(ii, 1), 500 * elems(ii, 1), 'NiTi_mono', ...
                q_cv_1(1) + rand*sigma_o, q_cv_1(2) + rand*sigma_o, q_cv_1(3) + rand*sigma_o, q_cv_1(4) + rand*sigma_o, ...
                1.0 + rand*sigma_e, 1.0 + rand*sigma_e, 1.0 + rand*sigma_e, ...
                rand*sigma_e, rand*sigma_e, rand*sigma_e);
    elseif(elems_phase(ii) == 3)
        fprintf(f, '%14.8f, %14.8f, %14.8f, %s, %14.8f, %14.8f, %14.8f, %14.8f, %14.8f, %14.8f, %14.8f, %14.8f, %14.8f, %14.8f\n', ...
                500 * elems(ii, 1), 500 * elems(ii, 1), 500 * elems(ii, 1), 'NiTi_mono', ...
                q_cv_2(1) + rand*sigma_o, q_cv_2(2) + rand*sigma_o, q_cv_2(3) + rand*sigma_o, q_cv_2(4) + rand*sigma_o, ...
                1.0 + rand*sigma_e, 1.0 + rand*sigma_e, 1.0 + rand*sigma_e, ...
                rand*sigma_e, rand*sigma_e, rand*sigma_e);
    end
end
fclose(f);
