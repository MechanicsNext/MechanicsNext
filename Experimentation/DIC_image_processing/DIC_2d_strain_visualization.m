% 2D DIC strain analysis
%
% 1. Run handles_ncorr = ncorr
% 2. Save the data and load the MAT file
% 3. Run this
%
num_of_dic_images = numel(current_save);
num_of_dic_structs = numel(data_dic_save);
num_of_dic_images_in_structs = zeros(num_of_dic_structs, 1);
for i = 1:num_of_dic_structs
    num_of_dic_images_in_structs(i) = numel(data_dic_save(i).strains);
end
%
% exx = 2D distribution of xx strain component in the ROI
exx = cell(num_of_dic_images,1);
% Save strain from each DIC result
dic_image_counter = 1;
for i = 1:num_of_dic_structs
    for j = 1:num_of_dic_images_in_structs(i)
        exx{dic_image_counter} = data_dic_save(i).strains(j).plot_exx_cur_formatted;
        dic_image_counter = dic_image_counter + 1;
    end
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
h0 = figure;
plot(exx_avg)
errorbar(exx_avg, exx_std)
xlabel('DIC frame')
ylabel('Strain (xx)')
%ylim([0 3e-3])
%
% Save an image sequence from 2d strain plots
make_image_sequence = 1;
strain_limits = [0.0 0.12];
plot_save_dir = '/Users/Harshad/Documents/Work/Research_data/HEDM_Study_of_Deformation_in_SMA/CHESS_2015_June_Steel/exx_2d';
padding = 0.4; % Padding around the plot. 0 < padding < 0.2
start_image_num = 1;
end_image_num = 845;
num_of_dic_images = 845;

if(make_image_sequence)
    close all;
    if(end_image_num > num_of_dic_images || end_image_num < 1)
        end_image_num = num_of_dic_images;
    end
    %
    for i = start_image_num:end_image_num
        % Get bounding box for the first frame. We will crop all the images
        % to that bbox
        if(i == 1 || i == start_image_num)
            exx_temp = exx{i};
            exx_temp = abs(floor(exx_temp/max(abs(exx_temp(:)))*256));
            exx_temp = uint8(exx_temp);
            % Get bounding box
            exx_bbox = regionprops(exx_temp ~= 0, 'BoundingBox');
            exx_bbox = exx_bbox.BoundingBox;
            % We dilate the bbox by 20% so that the deformed configuration
            % fits in it as well.
            exx_bbox(1) = floor(exx_bbox(1) - padding*exx_bbox(3));
            exx_bbox(2) = floor(exx_bbox(2) - padding*exx_bbox(4));
            exx_bbox(3) = floor(exx_bbox(3) + 2.0*padding*exx_bbox(3));
            exx_bbox(4) = floor(exx_bbox(4) + 2.0*padding*exx_bbox(4));
        end
        % Plotting
        exx_temp = exx{i};
        h = figure;
        % Plot per cent strain
        pcolor(rot90(exx_temp(exx_bbox(2):(exx_bbox(2) + exx_bbox(4)), ...
                     exx_bbox(1):(exx_bbox(1) + exx_bbox(3)))) * 100.0);
        shading flat;
        axis equal;
        axis off;
        colormap jet;
        colorbar;
        set(gca, 'FontSize', 16)
        caxis(100*[strain_limits(1) strain_limits(2)]);
        title(['Strain = ' sprintf('%6.2f', 100*exx_avg(i)) ' +- ' sprintf('%8.4f', 100*exx_std(i)) ' %']);
        print(h, '-dpng', fullfile(plot_save_dir, ['strain_2d_' sprintf('%04d', i)]));
        close(h);
    end
end