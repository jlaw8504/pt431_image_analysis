function [meta_lengths, ana_lengths] = plasmid_stretch_lengths_batch(threshold,folder)
%% Read in the Data
%select the folder
cd(folder);
%get the list of mat files
files = dir('*.mat');
%loop to get varibles into matrices
data_mat = [];
for n = 1:length(files)
    data{n} = load(files(n).name);
    data_mat = [data_mat; data{n}.data_cell(2:end,1:11)];
end

%% Seperate out the split data from the rest
counter = 1;
for i = 1:length(data_mat)
    if data_mat{i,11} == 0
        data_mat_no_split(counter,:) = data_mat(i,:);
        counter = counter + 1;
    end
end

%% Separate streched from foci data
stretch_count = 1;
foci_count = 1;
for j = 1:length(data_mat_no_split)
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
for k = 1:length(foci_data)
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
for k = 1:length(stretch_data)
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

%% Gather the distributions of stretch lengths
% >500 nm
if exist('stretch_less_500') == 1
    stretches_500 = cell2mat(stretch_less_500(:,2));
else
    stretches_500 = [];
end
% 500 - 1000 nm
if exist('stretch_less_1000') == 1
    stretches_1000 = cell2mat(stretch_less_1000(:,2));
else
    stretches_1000 = [];
end
% 1000 - 1500 nm
if exist('stretch_less_1500') == 1
    stretches_1500 = cell2mat(stretch_less_1500(:,2));
else
    stretches_1500 = [];
end
% 1500 - 2000 nm
if exist('stretch_less_2000') == 1
    stretches_2000 = cell2mat(stretch_less_2000(:,2));
else
    stretches_2000 = [];
end
% 2000 - 2500 nm
if exist('stretch_less_2500') == 1
    stretches_2500 = cell2mat(stretch_less_2500(:,2));
else
    stretches_2500 = [];
end
% 2500 - 3000
if exist('stretch_less_3000') == 1
    stretches_3000 = cell2mat(stretch_less_3000(:,2));
else
    stretches_3000 = [];
end
% 3000 - 3500
if exist('stretch_less_3500') == 1
    stretches_3500 = cell2mat(stretch_less_3500(:,2));
else
    stretches_3500 = [];
end
% 3500 - 4000
if exist('stretch_less_4000') == 1
    stretches_4000 = cell2mat(stretch_less_4000(:,2));
else
    stretches_4000 = [];
end
% >4000
if exist('stretch_greater_4000') == 1
    stretches_greater_4000 = cell2mat(stretch_greater_4000(:,2));
else
    stretches_greater_4000 = [];
end

%% Gather the metaphase lengths and anaphase lengths
meta_lengths = [stretches_1500;stretches_2000];
ana_lengths = [stretches_2500;stretches_3000;stretches_3500;stretches_4000];
%% calculate the metaphase and anaphase stretch freq
%If variable doesn't exist then size is 0
if exist('stretch_less_1500','var')== 1
    num_stretch_1500 = size(stretch_less_1500);
else
    num_stretch_1500 = 0;
end
if exist('stretch_less_2000','var')== 1
    num_stretch_2000 = size(stretch_less_2000);
else
    num_stretch_2000 = 0;
end
if exist('stretch_less_2500','var')== 1
    num_stretch_2500 = size(stretch_less_2500);
else
    num_stretch_2500 = 0;
end
if exist('stretch_less_3000','var')== 1
    num_stretch_3000 = size(stretch_less_3000);
else
    num_stretch_3000 = 0;
end
if exist('stretch_less_3500','var')== 1
    num_stretch_3500 = size(stretch_less_3500);
else
    num_stretch_3500 = 0;
end
if exist('stretch_less_4000','var')== 1
    num_stretch_4000 = size(stretch_less_4000);
else
    num_stretch_4000 = 0;
end
if exist('foci_less_1500','var')== 1
    num_foci_1500 = size(foci_less_1500);
else
    num_foci_1500 = 0;
end
if exist('foci_less_2000','var')== 1
    num_foci_2000 = size(foci_less_2000);
else
    num_foci_2000 = 0;
end
if exist('foci_less_2500','var')== 1
    num_foci_2500 = size(foci_less_2500);
else
    num_foci_2500 = 0;
end
if exist('foci_less_3000','var')== 1
    num_foci_3000 = size(foci_less_3000);
else
    num_foci_3000 = 0;
end
if exist('foci_less_3500','var')== 1
    num_foci_3500 = size(foci_less_3500);
else
    num_foci_3500 = 0;
end
if exist('foci_less_4000','var')== 1
    num_foci_4000 = size(foci_less_4000);
else
    num_foci_4000 = 0;
end
freq_meta = (num_stretch_1500 + num_stretch_2000)/...
    (num_stretch_1500 + num_foci_1500 + num_stretch_2000 + ...
    num_foci_2000);
freq_ana = (num_stretch_2500 + num_stretch_3000 + num_stretch_3500 + ...
    num_stretch_4000) / (num_stretch_2500 + num_foci_2500 + num_stretch_3000 + ...
    num_foci_3000 + num_stretch_3500 + num_foci_3500 + num_stretch_4000 + ...
    num_foci_4000);
%% Show the means and stds
freq_meta
meta_mean = mean(meta_lengths)
meta_std = std(meta_lengths)
freq_ana
ana_mean = mean(ana_lengths)
ana_std = std(ana_lengths)
end