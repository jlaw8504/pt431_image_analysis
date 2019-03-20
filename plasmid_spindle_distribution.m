function [spindles_total, spindles_mean, spindles_std] = plasmid_spindle_distribution
%This function will calculate the correlation coefficient of the major axis
%length of stretched plasmids (based on provided aspect ratio threshold)
%to 3D spindle length
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

%% Calculate and plot histograms
% gather spindle data for all data
spindles_total = cell2mat(data_mat_no_split(:,6));
spindles_mean = mean(spindles_total);
spindles_std = std(spindles_total);

