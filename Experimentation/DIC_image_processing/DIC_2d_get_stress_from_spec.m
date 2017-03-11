spec_file_dir = '/Users/Harshad/Documents/Work/Research_data/HEA_HEDM_CHESS_2015/';
spec_log_file_name = 'spec.log';
cross_section_area = 1;

log_file = fopen(fullfile(spec_file_dir, spec_log_file_name), 'r');
image_file_numbers = [];
force = [];
displacement = [];
newimages = [];
ff_image_points = [];
start_scan_flag = 1;
ff_flag = 0;

% Extract load/displacement from log.
while(~feof(log_file))
    log_line_next = fgetl(log_file);
    if(~isempty(findstr(log_line_next, 'scan_ome_ff')) || ~isempty(findstr(log_line_next, 'slew_ome')))
        ff_flag = 1;
    end
    %if(start_scan_flag == 0 && ~isempty(findstr(log_line, '#D Thu Oct 30')))
    %    start_scan_flag = 1;
    %end
    if(~isempty(findstr(log_line_next, 'force=')) && ...
            ~isempty(findstr(log_line_next, 'displacement=')))
        log_line = fgetl(log_file);
        if(start_scan_flag ~= 0 && ~isempty(findstr(log_line, 'DIC image saved')))
            for ii = 13:-1:7
                if isempty(str2num(log_line(end - ii))) == 1
                    imindex = ii-1;
                end
            end
            image_file_numbers(end+1) = str2num(log_line(end-imindex:end-6));
            newimages(end+1,1) = 1e5 + image_file_numbers(end);
            split_log_line = strsplit(log_line_next, '=');
            split_log_line_2 = strsplit(mat2str(cell2mat(split_log_line(2))), ',');
            force_str_temp = mat2str(cell2mat(split_log_line_2(1)));
            force(end+1) = str2num(force_str_temp(4:end-1));
            split_log_line_3 = mat2str(cell2mat(split_log_line(3)));
            split_log_line_4 = split_log_line_3(2:end-2);
            %displacement_str_temp = mat2str(cell2mat(split_log_line_4));
            displacement(end+1) = str2num(split_log_line_4);
            %
            if(ff_flag == 0)
                ff_image_points(end+1) = 0;
            else
                ff_image_points(end+1) = 1;
                ff_flag = 0;
            end
        end
    end
end
fclose(log_file);
% Plot for fun
plot(displacement-displacement(1), force)
xlabel('Relative crosshead displacement (mm)');
ylabel('Force (N)');
%
% dic_files = dir(DIC_image_dir);
% dic_file_names = {dic_files.name};
% force_filtered = zeros(numel(dic_file_names) - 2, 1);
% displacement_filtered = zeros(numel(dic_file_names) - 2, 1);
% image_file_numbers_filtered = zeros(numel(dic_file_names) - 2, 1);
% for i = 3:(numel(dic_file_names) - 2)
%     if(i > numel(force))
%         exit;
%     end
%     fname = dic_file_names{i};
%     fname_2 = strsplit(fname, '_');
%     fname_3 = strsplit(fname_2{2}, '.');
%     fnumber = str2num(fname_3{1});
%     disp(fnumber)
%     force_filtered(i - 2) = force(fnumber);
%     displacement_filtered(i - 2) = displacement(fnumber);
%     image_file_numbers_filtered(i - 2) = image_file_numbers(fnumber);
% end