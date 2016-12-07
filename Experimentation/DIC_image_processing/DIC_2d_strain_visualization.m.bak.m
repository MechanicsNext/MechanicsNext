% 2D DIC strain analysis
%
% 1. Run handles_ncorr = ncorr
% 2. Save the data and load the MAT file
% 3. Run this
%
num_of_dic_images = numel(current_save);
%
% exx = 2D distribution of xx strain component in the ROI
exx = cell(num_of_dic_images,1);
% Save strain from each DIC result
for i = 1:num_of_dic_images
    exx{i} = data_dic_save.strains(i).plot_exx_cur_formatted;
end
% Calculate mean and std deviation for each frame
% Calculate for exx now
exx_avg = zeros(num_of_dic_images,1);
exx_std = zeros(num_of_dic_images,1);
for i = 1:num_of_dic_images
    exx_temp = exx{i};
    exx_avg(i) = mean(exx_temp(exx_temp ~= 0));
    exx_std(i) = std(exx_temp(exx_temp ~= 0));
end
% Plot avg strain and errorbars
plot(exx_avg)
errorbar(exx_avg, exx_std)
xlabel('DIC frame')
ylabel('Strain (xx)')
%ylim([0 3e-3])
