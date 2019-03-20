function [spindles_foci, spindles_foci_mean, spindles_foci_std]...
    = plasmid_foci_spindle_distribution(threshold)
%This function will return the spindle lengths of cells with stretched GFP
%signals
%% Read in the Data
%select the folder
folder = uigetdir;
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
%% Calculate and plot histograms
% gather spindle data for all data
spindles_foci = cell2mat(foci_data(:,6));
spindles_foci_mean = mean(spindles_foci);
spindles_foci_std = std(spindles_foci);

