% Where are ncorr result files located
strain_result_mat_directory = '/Users/Harshad/Documents/Work/Research_data/HEA_HEDM_CHESS_2015/DIC_data/';
% List of ncorr result mat files in the correct order
strain_result_mat_files = { ...
    'strain_2d_1', ...
    'strain_2d_2', ...
    };
ignore_these_many_DIC_results_at_the_start = 0;
ignore_these_many_DIC_results_at_the_end = 0;
% Aggregate the strain stuff
exx_avg_aggregate = [];
exx_std_aggregate = [];
for i = 1:numel(strain_result_mat_files)
    mat_file_name_i = fullfile(strain_result_mat_directory, [strain_result_mat_files{i} '.mat']);
    load(mat_file_name_i);
    run DIC_2d_strain_visualization;
    exx_avg_aggregate = [exx_avg_aggregate; exx_avg];
    exx_std_aggregate = [exx_std_aggregate; exx_std];
end
%
exx_avg_aggregate = exx_avg_aggregate((ignore_these_many_DIC_results_at_the_start+1):end);
exx_avg_aggregate = exx_avg_aggregate(1:(end-ignore_these_many_DIC_results_at_the_end));
% Get forces from spec
% Make sure that the DIC_2d_get_stress_from_spec.m script has 
% correct spec location info
% Make sure ncorr results and spec log ends at the same stage
% THIS WHOLE THING ONLY WORKS IF ONLY THE DIC FILES AT THE 
% BEGINNING OF THE SERIES ARE JUNK. IF THERE IS JUNK DATA
% IN BETWEEN THIS FAILS.
run DIC_2d_get_stress_from_spec;
force_filtered = force((end - numel(exx_avg_aggregate) + 1):end);
ff_image_points_filtered = ff_image_points((end - numel(exx_avg_aggregate) + 1):end);
% Make plot
plot(exx_avg_aggregate', force_filtered);
xlabel('Strain (X)');
ylabel('Load');
