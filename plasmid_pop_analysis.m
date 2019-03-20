function [ table ] = plasmid_pop_analysis( threshold )

%% Read in the Data
%select the folder
folder = uigetdir;
cd(folder);
%get the list of mat files
files = dir('*.mat');
%loop to get varibles into matrices
data_mat = [];
for n = 1:size(files,1)
    data{n} = load(files(n).name);
    data_mat = [data_mat; data{n}.data_cell(2:end,1:11)];
end

%% Seperate out the split data from the rest
counter = 1;
for i = 1:size(data_mat,1)
    if data_mat{i,11} == 0
        data_mat_no_split(counter,:) = data_mat(i,:);
        counter = counter + 1;
    end
end

%% Separate streched from foci data
stretch_count = 1;
foci_count = 1;
for j = 1:size(data_mat_no_split,1)
    if data_mat_no_split{j,4} > threshold
        stretch_data(stretch_count,:) = data_mat_no_split(j,:);
        stretch_count = stretch_count + 1;
    else
        foci_data(foci_count,:) = data_mat_no_split(j,:);
        foci_count = foci_count + 1;
    end
end

%% Separate the foci into bins based on spindle length
%set up counters
counter_500 = 1;
counter_1000 = 1;
counter_1500 = 1;
counter_2000 = 1;
counter_2500 = 1;
counter_3000 = 1;
counter_3500 = 1;
counter_4000 = 1;
counter_greater = 1;
%loop through foci data and assign to bins
for k = 1:size(foci_data,1)
    if foci_data{k,6} < 500
        foci_less_500(counter_500,:) = foci_data(k,:);
        counter_500 = counter_500 + 1;
    elseif foci_data{k,6} < 1000
        foci_less_1000(counter_1000,:) = foci_data(k,:);
        counter_1000 = counter_1000 + 1;
    elseif foci_data{k,6} < 1500
        foci_less_1500(counter_1500,:) = foci_data(k,:);
        counter_1500 = counter_1500 + 1;
    elseif foci_data{k,6} < 2000
        foci_less_2000(counter_2000,:) = foci_data(k,:);
        counter_2000 = counter_2000 + 1;
    elseif foci_data{k,6} < 2500
        foci_less_2500(counter_2500,:) = foci_data(k,:);
        counter_2500 = counter_2500 + 1;
    elseif foci_data{k,6} < 3000
        foci_less_3000(counter_3000,:) = foci_data(k,:);
        counter_3000 = counter_3000 + 1;
    elseif foci_data{k,6} < 3500
        foci_less_3500(counter_3500,:) = foci_data(k,:);
        counter_3500 = counter_3500 + 1;
    elseif foci_data{k,6} < 4000
        foci_less_4000(counter_4000,:) = foci_data(k,:);
        counter_4000 = counter_4000 + 1;
    else
        foci_greater_4000{counter_greater,:} = foci_data(k,:);
        counter_greater = counter_greater + 1;
    end
end

%% Separate the stretched signals into bins based on spindle length
%reset the counters
counter_500 = 1;
counter_1000 = 1;
counter_1500 = 1;
counter_2000 = 1;
counter_2500 = 1;
counter_3000 = 1;
counter_3500 = 1;
counter_4000 = 1;
counter_greater = 1;
%loop through stretch data and assign to bins
for k = 1:size(stretch_data,1)
    if stretch_data{k,6} < 500
        stretch_less_500(counter_500,:) = stretch_data(k,:);
        counter_500 = counter_500 + 1;
    elseif stretch_data{k,6} < 1000
        stretch_less_1000(counter_1000,:) = stretch_data(k,:);
        counter_1000 = counter_1000 + 1;
    elseif stretch_data{k,6} < 1500
        stretch_less_1500(counter_1500,:) = stretch_data(k,:);
        counter_1500 = counter_1500 + 1;
    elseif stretch_data{k,6} < 2000
        stretch_less_2000(counter_2000,:) = stretch_data(k,:);
        counter_2000 = counter_2000 + 1;
    elseif stretch_data{k,6} < 2500
        stretch_less_2500(counter_2500,:) = stretch_data(k,:);
        counter_2500 = counter_2500 + 1;
    elseif stretch_data{k,6} < 3000
        stretch_less_3000(counter_3000,:) = stretch_data(k,:);
        counter_3000 = counter_3000 + 1;
    elseif stretch_data{k,6} < 3500
        stretch_less_3500(counter_3500,:) = stretch_data(k,:);
        counter_3500 = counter_3500 + 1;
    elseif stretch_data{k,6} < 4000
        stretch_less_4000(counter_4000,:) = stretch_data(k,:);
        counter_4000 = counter_4000 + 1;
    else
        stretch_greater_4000(counter_greater,:) = stretch_data(k,:);
        counter_greater = counter_greater + 1;
    end
end

%% Calculate the stretch frequency by spindle size
% setup the number of stretch and foci per spindle length
% 500 nm or less
if exist('stretch_less_500','var') == 1 && exist('foci_less_500','var') == 1
    freq_500 = size(stretch_less_500,1)/(size(stretch_less_500,1) + ...
        size(foci_less_500,1));
elseif exist('stretch_less_500','var') == 0 && exist('foci_less_500','var') == 1
    freq_500 = 0;
elseif exist('stretch_less_500','var') == 1 && exist('foci_less_500','var') == 0
    freq_500 = 1;
else
    freq_500 = [];
end
% 500-1000 nm
if exist('stretch_less_1000','var') == 1 && exist('foci_less_1000','var') == 1
    freq_1000 = size(stretch_less_1000,1)/(size(stretch_less_1000,1) + ...
        size(foci_less_1000,1));
elseif exist('stretch_less_1000','var') == 0 && exist('foci_less_1000','var') == 1
    freq_1000 = 0;
elseif exist('stretch_less_1000','var') == 1 && exist('foci_less_1000','var') == 0
    freq_1000 = 1;
else
    freq_1000 = [];
end
% 1000-1500 nm
if exist('stretch_less_1500','var') == 1 && exist('foci_less_1500','var') == 1
    freq_1500 = size(stretch_less_1500,1)/(size(stretch_less_1500,1) + ...
        size(foci_less_1500,1));
elseif exist('stretch_less_1500','var') == 0 && exist('foci_less_1500','var') == 1
    freq_1500 = 0;
elseif exist('stretch_less_1500','var') == 1 && exist('foci_less_1500','var') == 0
    freq_1500 = 1;
else
    freq_1500 = [];
end
% 1500-2000 nm
if exist('stretch_less_2000','var') == 1 && exist('foci_less_2000','var') == 1
    freq_2000 = size(stretch_less_2000,1)/(size(stretch_less_2000,1) + ...
        size(foci_less_2000,1));
elseif exist('stretch_less_2000','var') == 0 && exist('foci_less_2000','var') == 1
    freq_2000 = 0;
elseif exist('stretch_less_2000','var') == 1 && exist('foci_less_2000','var') == 0
    freq_2000 = 1;
else
    freq_2000 = [];
end
% 2000-2500 nm
if exist('stretch_less_2500','var') == 1 && exist('foci_less_2500','var') == 1
    freq_2500 = size(stretch_less_2500,1)/(size(stretch_less_2500,1) + ...
        size(foci_less_2500,1));
elseif exist('stretch_less_2500','var') == 0 && exist('foci_less_2500','var') == 1
    freq_2500 = 0;
elseif exist('stretch_less_2500','var') == 1 && exist('foci_less_2500','var') == 0
    freq_2500 = 1;
else
    freq_2500 = [];
end
% 2500-3000 nm
if exist('stretch_less_3000','var') == 1 && exist('foci_less_3000','var') == 1
    freq_3000 = size(stretch_less_3000,1)/(size(stretch_less_3000,1) + ...
        size(foci_less_3000,1));
elseif exist('stretch_less_3000','var') == 0 && exist('foci_less_3000','var') == 1
    freq_3000 = 0;
elseif exist('stretch_less_3000','var') == 1 && exist('foci_less_3000','var') == 0
    freq_3000 = 1;
else
    freq_3000 = [];
end
% 3000-3500 nm
if exist('stretch_less_3500','var') == 1 && exist('foci_less_3500','var') == 1
    freq_3500 = size(stretch_less_3500,1)/(size(stretch_less_3500,1) + ...
        size(foci_less_3500,1));
elseif exist('stretch_less_3500','var') == 0 && exist('foci_less_3500','var') == 1
    freq_3500 = 0;
elseif exist('stretch_less_3500','var') == 1 && exist('foci_less_3500','var') == 0
    freq_3500 = 1;
else
    freq_3500 = [];
end
% 3500-4000 nm
if exist('stretch_less_4000','var') == 1 && exist('foci_less_4000','var') == 1
    freq_4000 = size(stretch_less_4000,1)/(size(stretch_less_4000,1) + ...
        size(foci_less_4000,1));
elseif exist('stretch_less_4000','var') == 0 && exist('foci_less_4000','var') == 1
    freq_4000 = 0;
elseif exist('stretch_less_4000','var') == 1 && exist('foci_less_4000','var') == 0
    freq_4000 = 1;
else
    freq_4000 = [];
end
% >4000 nm
if exist('stretch_greater_4000','var') == 1 && exist('foci_greater_4000','var') == 1
    freq_greater_4000 = size(stretch_greater_4000,1)/(size(stretch_greater_4000,1) + ...
        size(foci_greater_4000,1));
elseif exist('stretch_greater_4000','var') == 0 && exist('foci_greater_4000','var') == 1
    freq_greater_4000 = 0;
elseif exist('stretch_greater_4000','var') == 1 && exist('foci_greater_4000','var') == 0
    freq_greater_4000 = 1;
else
    freq_greater_4000 = [];
end

%% Calculate the average plasmid stretch with standard deviation
% >500 nm
if exist('stretch_less_500') == 1
    mean_500 = mean(cell2mat(stretch_less_500(:,2)));
    std_500 = std(cell2mat(stretch_less_500(:,2)));
else
    mean_500 = [];
    std_500 = [];
end
% 500 - 1000 nm
if exist('stretch_less_1000') == 1
    mean_1000 = mean(cell2mat(stretch_less_1000(:,2)));
    std_1000 = std(cell2mat(stretch_less_1000(:,2)));
else
    mean_1000 = [];
    std_1000 = [];
end
% 1000 - 1500 nm
if exist('stretch_less_1500') == 1
    mean_1500 = mean(cell2mat(stretch_less_1500(:,2)));
    std_1500 = std(cell2mat(stretch_less_1500(:,2)));
else
    mean_1500 = [];
    std_1500 = [];
end
% 1500 - 2000 nm
if exist('stretch_less_2000') == 1
    mean_2000 = mean(cell2mat(stretch_less_2000(:,2)));
    std_2000 = std(cell2mat(stretch_less_2000(:,2)));
else
    mean_2000 = [];
    std_2000 = [];
end
% 2000 - 2500 nm
if exist('stretch_less_2500') == 1
    mean_2500 = mean(cell2mat(stretch_less_2500(:,2)));
    std_2500 = std(cell2mat(stretch_less_2500(:,2)));
else
    mean_2500 = [];
    std_2500 = [];
end
% 2500 - 3000
if exist('stretch_less_3000') == 1
    mean_3000 = mean(cell2mat(stretch_less_3000(:,2)));
    std_3000 = std(cell2mat(stretch_less_3000(:,2)));
else
    mean_3000 = [];
    std_3000 = [];
end
% 3000 - 3500
if exist('stretch_less_3500') == 1
    mean_3500 = mean(cell2mat(stretch_less_3500(:,2)));
    std_3500 = std(cell2mat(stretch_less_3500(:,2)));
else
    mean_3500 = [];
    std_3500 = [];
end
% 3500 - 4000
if exist('stretch_less_4000') == 1
    mean_4000 = mean(cell2mat(stretch_less_4000(:,2)));
    std_4000 = std(cell2mat(stretch_less_4000(:,2)));
else
    mean_4000 = [];
    std_4000 = [];
end
% >4000
if exist('stretch_greater_4000') == 1
    mean_greater_4000 = mean(cell2mat(stretch_greater_4000(:,2)));
    std_greater_4000 = std(cell2mat(stretch_greater_4000(:,2)));
else
    mean_greater_4000 = [];
    std_greater_4000 = [];
end

%% Calculate the number of plasmids in each bin
% <500 nm
if exist('stretch_less_500','var') == 1
    num_stretch_500 = size(stretch_less_500,1);
else
    num_stretch_500 = 0;
end
if exist('foci_less_500','var') == 1
    num_foci_500 = size(foci_less_500,1);
else
    num_foci_500 = 0;
end
num_total_500 = num_stretch_500 + num_foci_500;

% 500 - 1000 nm
if exist('stretch_less_1000','var') == 1
    num_stretch_1000 = size(stretch_less_1000,1);
else
    num_stretch_1000 = 0;
end
if exist('foci_less_1000','var') == 1
    num_foci_1000 = size(foci_less_1000,1);
else
    num_foci_1000 = 0;
end
num_total_1000 = num_stretch_1000 + num_foci_1000;

% 1000 - 1500 nm
if exist('stretch_less_1500','var') == 1
    num_stretch_1500 = size(stretch_less_1500,1);
else
    num_stretch_1500 = 0;
end
if exist('foci_less_1500','var') == 1
    num_foci_1500 = size(foci_less_1500,1);
else
    num_foci_1500 = 0;
end
num_total_1500 = num_stretch_1500 + num_foci_1500;

% 1500 - 2000 nm
if exist('stretch_less_2000','var') == 1
    num_stretch_2000 = size(stretch_less_2000,1);
else
    num_stretch_2000 = 0;
end
if exist('foci_less_2000','var') == 1
    num_foci_2000 = size(foci_less_2000,1);
else
    num_foci_2000 = 0;
end
num_total_2000 = num_stretch_2000 + num_foci_2000;

% 2000 - 2500 nm
if exist('stretch_less_2500','var') == 1
    num_stretch_2500 = size(stretch_less_2500,1);
else
    num_stretch_2500 = 0;
end
if exist('foci_less_2500','var') == 1
    num_foci_2500 = size(foci_less_2500,1);
else
    num_foci_2500 = 0;
end
num_total_2500 = num_stretch_2500 + num_foci_2500;

% 2000 - 2500 nm
if exist('stretch_less_3000','var') == 1
    num_stretch_3000 = size(stretch_less_3000,1);
else
    num_stretch_3000 = 0;
end
if exist('foci_less_3000','var') == 1
    num_foci_3000 = size(foci_less_3000,1);
else
    num_foci_3000 = 0;
end
num_total_3000 = num_stretch_3000 + num_foci_3000;

% 3000 - 3500 nm
if exist('stretch_less_3500','var') == 1
    num_stretch_3500 = size(stretch_less_3500,1);
else
    num_stretch_3500 = 0;
end
if exist('foci_less_3500','var') == 1
    num_foci_3500 = size(foci_less_3500,1);
else
    num_foci_3500 = 0;
end
num_total_3500 = num_stretch_3500 + num_foci_3500;

% 3500 - 4000 nm
if exist('stretch_less_4000','var') == 1
    num_stretch_4000 = size(stretch_less_4000,1);
else
    num_stretch_4000 = 0;
end
if exist('foci_less_4000','var') == 1
    num_foci_4000 = size(foci_less_4000,1);
else
    num_foci_4000 = 0;
end
num_total_4000 = num_stretch_4000 + num_foci_4000;

% > 4000 nm
if exist('stretch_greater_4000','var') == 1
    num_stretch_greater_4000 = size(stretch_greater_4000,1);
else
    num_stretch_greater_4000 = 0;
end
if exist('foci_greater_4000','var') == 1
    num_foci_greater_4000 = size(foci_greater_4000,1);
else
    num_foci_greater_4000 = 0;
end
num_total_greater_4000 = num_stretch_greater_4000 + num_foci_greater_4000;

%% Summary stretch freqency and magnitude in a table
spindle_bins = {'<500'; '500 - 1000'; '1000 - 1500'; '1500 - 2000'; ...
    '2000 - 2500'; '2500 - 3000'; '3000 - 3500'; '3500 - 4000'; '>4000'};
numbers = {num_total_500, num_total_1000, num_total_1500,...
    num_total_2000, num_total_2500, num_total_3000, num_total_3500, ...
    num_total_4000, num_total_greater_4000}';
freqs = {freq_500, freq_1000, freq_1500, freq_2000, freq_2500, ...
    freq_3000, freq_3500, freq_4000, freq_greater_4000}';
means = {mean_500, mean_1000, mean_1500, mean_2000, mean_2500, ...
    mean_3000, mean_3500, mean_4000, mean_greater_4000}';
stds = {std_500, std_1000, std_1500, std_2000, std_2500, ...
    std_3000, std_3500, std_4000, std_greater_4000}';
summary_cell = [spindle_bins, numbers, freqs, means, stds];
table = cell2table(summary_cell,'VariableNames',{'SpindleLength'...
    'N' 'StretchFreqency' 'MeanStretchLength'...
    'StdDevStretchLength'});
end