load image_list.mat

source_dir = '../../svt/svt1/img/';
results_dir = 'results_second_run/';

if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end
% for i=1:size(image_list, 1)
for i=1:100
    source_file = [source_dir image_list{i}];
    result_file = [results_dir 'result' image_list{i}];
    result = detecttext(source_file);
    imwrite(result, result_file, 'jpg');
    
end
