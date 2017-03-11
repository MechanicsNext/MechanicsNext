function [] = print_rotated_twin_habit_plane_normals(TEXTURE, loading_axis, foil_plane, print_flag, print_sort_order, plot_flag, U)
% SMA Calculations: Habit/Twin Elements in the global coordinates
%
% Rotate twin elements (a, n) and habit plane elements (b, m) from crystal
% to global coordinate system.
% DEPENDENCIES: MTEX toolbox, crystallography routines by Harshad Paranjape
% INPUTS:
% TEXTURE      = crystal orientation of the grain. Quaternion. See heXRD docs
%                for the explanation of quaternion.
% loading_axis = orientation along which load is applied. In global
%                coordinates
% foil_plane   = the surface normal for the sample. This is the plane on
%                which m, n are projected for plotting. DIC images are taken along the
% foil_plane   = Twin/habit plane projections will be done on this plane
% print_flag   = 1 prints the a, n, b, m values in the crystal coordinate
%                system. 0 does not print anything
% Sort order   = Order in which printed b, m values are sorted. 
%                0 | 1: 0 = sort by CV number, 1 = descending by tstrain
% plot_flag    = 1 plots 2D projection of m, n. Look at the code at the bottom
%                to choose specific h.p.v.s for plotting. You don't want this to make 192
%                plots at once for NiTi.
% U            = Cell array of transformation deformation gradients (not strain)
% Get a, n, b, m
[~, ~, twin_elements, am_elements, corresponding_twin] = habit_calculation_3d(U);
% Rotate everything from crystal to global coordinates
for itex=1:size(TEXTURE, 1)
    % euler = TEXTURE(itex, :)*pi/180.0;
    % Get orientation matrix (crystal to sample!).
    % tlg=euler2g(euler);
    tlg = matrix(rotation(quaternion(TEXTURE(itex,:)')));
    % Rotate b, m for habit planes
    cb0 = am_elements(:, 3:5)';
    cm0 = am_elements(:, 6:8)';
    % Transform to sample
    b0 = tlg*cb0;
    m0 = tlg*cm0;
    theta_habit = zeros(size(m0, 2), 1);
    % Angle between projection of habit plane **normal** on foil plane and
    % the loading axis
    for i=1:size(theta_habit, 1)
        m0_projected = m0(:, i)' - foil_plane*(dot(m0(:, i)', foil_plane)/dot(foil_plane, foil_plane));
        theta_habit(i) = acosd(dot(m0_projected, loading_axis)/sqrt(norm(m0_projected)*norm(loading_axis)));
        if(theta_habit(i) > 90)
            theta_habit(i) = theta_habit(i) - 180.0;
        end
        % Add 90 degrees to get the angle made by habit plane **trace**
        % with the loading axis
        theta_habit(i) = theta_habit(i) + 90;
    end
    schmid_habit = zeros(size(m0, 2), 1);
    for i=1:size(theta_habit, 1)
        schmid_habit(i) = loading_axis*(0.5*(kron(b0(:,i), m0(:, i)')+kron(m0(:,i), b0(:, i)')))*loading_axis';
    end
    am_elements_global = [am_elements(:, 1:2) b0' m0' am_elements(:, 9)];
    % Rotate a, n for twin planes
    ca0 = twin_elements(:, 3:5)';
    cn0 = twin_elements(:, 6:8)';
    % Transform to sample
    a0 = tlg*ca0;
    n0 = tlg*cn0;
    theta_twin = zeros(size(n0, 2), 1);
    % Angle between projection of twin plane **normal** on foil plane and
    % the loading axis
    for i=1:size(theta_twin, 1)
        n0_projected = n0(:, i)' - foil_plane*(dot(n0(:, i)', foil_plane)/dot(foil_plane, foil_plane));
        theta_twin(i) = acosd(dot(n0_projected, loading_axis)/sqrt(norm(n0_projected)*norm(loading_axis)));
        if(theta_twin(i) > 90)
            theta_twin(i) = theta_twin(i) - 180.0;
        end
        % Add 90 degrees to get the angle made by twin plane **trace** with
        % the loading axis
        theta_twin(i) = theta_twin(i) + 90;
    end
    twin_elements_global = [twin_elements(:, 1:2) a0' n0'];
    %% Print table
    if(print_flag == 1)
        % Print a, n, b, m
        % Print M-M twin elements
        twin_plane_data_global = [[1:size(twin_elements_global, 1)]' twin_elements_global theta_twin];
        invariant_plane_data_global = [[1:size(am_elements_global, 1)]' am_elements_global 100*schmid_habit theta_habit corresponding_twin'];
        if(print_sort_order == 0)
            disp('             CV1 CV2        |------------------a------------------|        |------------------n------------------|       Angle')
            disp(sprintf(' %4d\t %2d\t %2d\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %7.1f\n', [[1:size(twin_elements_global, 1)]' twin_elements_global theta_twin]'))
            % Print A-M invariant plane elements
            disp('       CV1  CV2     |------------------b------------------|        |------------------m------------------|         lambda                  Schmid          Angle')
            disp(sprintf(' %3d\t %2d\t %2d\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %7.1f\t %3.1f\t %2d\n', ...
                [[1:size(am_elements_global, 1)]' am_elements_global 100*schmid_habit theta_habit corresponding_twin']'))
        elseif(print_sort_order == 1)
            disp('Twin plane elements (a = direction, n = twin plane normal)')
            disp('         CV1     CV2        |------------------a------------------|        |------------------n------------------|       Angle')
            disp(sprintf(' %4d\t %2d\t %2d\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %7.1f\n', [[1:size(twin_elements_global, 1)]' twin_elements_global theta_twin]'))
            [~, idx] = sort(invariant_plane_data_global(:, 11));
            invariant_plane_data_global_sorted = invariant_plane_data_global(idx, :);
            % Print B-M invariant plane elements
            disp('Invariant (habit) plane elements (b = direction, m = habit plane normal)')
            disp('         CV1     CV2     |------------------b------------------|        |------------------m------------------|          lambda             Schmid      Angle')
            disp(sprintf(' %3d\t %2d\t %2d\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %7.1f\t %3.1f\t %2d\n', ...
                invariant_plane_data_global_sorted(:, :)'))
        end
    end
    % Plot AM and MM interface projection in foil_plane
    if(plot_flag == 1)
        am_elements_global(am_elements_global(:, 1) ~= 4, :) = [];
        am_elements_global(am_elements_global(:, 2) ~= 2, :) = [];
        twin_elements_global(twin_elements_global(:, 1) ~= 4, :) = [];
        twin_elements_global(twin_elements_global(:, 2) ~= 2, :) = [];
        %
        twin_trace_1 = cross(twin_elements_global(1, 6:8), foil_plane)/sqrt(norm(cross(twin_elements_global(1, 6:8), foil_plane)));
        twin_trace_2 = cross(twin_elements_global(2, 6:8), foil_plane)/sqrt(norm(cross(twin_elements_global(2, 6:8), foil_plane)));
        figure('Position', [200 200 1200 600])
        for i=1:8
            h = subaxis(2, 4, i, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05);
            habit_trace = cross(am_elements_global(i, 6:8), foil_plane)/sqrt(norm(cross(am_elements_global(i, 6:8), foil_plane)));
            hold on;
            line([-habit_trace(1) habit_trace(1)], [-habit_trace(2) habit_trace(2)], 'Color', 'g', 'LineWidth', 5);
            line([-twin_trace_1(1) twin_trace_1(1)], [-twin_trace_1(2) twin_trace_1(2)], 'Color', 'm', 'LineStyle', '-.', 'LineWidth', 5);
            line([-twin_trace_2(1) twin_trace_2(1)], [-twin_trace_2(2) twin_trace_2(2)], 'Color', 'm', 'LineStyle', ':', 'LineWidth', 5);
            for j=1:6
                slip_systems = [ 0  1  1  1  0  0; ...
                    0  1 -1  1  0  0; ...
                    1  0  1  0  1  0; ...
                    1  0 -1  0  1  0; ...
                    1  1  0  0  0  1; ...
                    -1  1  0  0  0  1];
                slip_trace = cross(slip_systems(j, 1:3), foil_plane)/sqrt(norm(cross(slip_systems(j, 1:3), foil_plane)));
                %line([-slip_trace(1) slip_trace(1)], [-slip_trace(2) slip_trace(2)], 'Color', rand(1, 3), 'LineStyle', ':', 'LineWidth', 2);
            end
            hold off;
            xlim([-1 1]); ylim([-1 1]);
            %        legend('Invariant plane', 'Type I twin', 'Type II twin', '0  1  1  1  0  0', '0  1 -1  1  0  0', '1  0  1  0  1  0', '1  0 -1  0  1  0', '1  1  0  0  0  1', '-1  1  0  0  0  1');
            legend('Invariant plane', 'Type I twin', 'Type II twin');
            axis square;
            set(h, 'XTickLabel', []);
            set(h, 'YTickLabel', []);
        end
    end
end
end
