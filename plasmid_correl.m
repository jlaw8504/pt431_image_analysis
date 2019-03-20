function [meta_stretch_corr_sq, ana_stretch_corr_sq] = plasmid_correl(threshold)
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

%% Separate streched from foci data
stretch_count = 1;
for j = 1:length(data_mat_no_split)
    if data_mat_no_split{j,4} > threshold
        stretch_data(stretch_count,1) = data_mat_no_split{j,2};
        stretch_data(stretch_count,2) = data_mat_no_split{j,6};
        stretch_count = stretch_count + 1;
    end
end
%% Using logical indexing to separate metaphase from anaphase spindles
meta_stretch = stretch_data(stretch_data(:,2) > 1000 ...
& stretch_data(:,2) < 2000,:);
ana_stretch = stretch_data(stretch_data(:,2) > 2000 ...
& stretch_data(:,2) < 4000,:);
%% calculate correlations^2 of stretched plasmids to spindle length
meta_stretch_corr_sq = corr(meta_stretch(:,1), meta_stretch(:,2))^2;
ana_stretch_corr_sq = corr(ana_stretch(:,1), ana_stretch(:,2))^2;
%% plot metaphase
figure;
plot(meta_stretch(:,1), meta_stretch(:,2), 'bx');
title('Metaphase Cells');
xlabel('Extended Plamsid Signal Length (nm)');
ylabel('3D Spindle Length (nm)');

%% calculate linear regression for metaphase and plot
slope = meta_stretch(:,1)\meta_stretch(:,2);
yCalc1 = slope*meta_stretch(:,1);
hold on;
plot(meta_stretch(:,1),yCalc1,'-k');
hold off;
%% plot anaphase
figure;
plot(ana_stretch(:,1), ana_stretch(:,2), 'bx');
title('Anaphase Cells');
xlabel('Extended Plamsid Signal Length (nm)');
ylabel('3D Spindle Length (nm)');

%% calculate linear regression for anaphase and plot
slope = ana_stretch(:,1)\ana_stretch(:,2);
yCalc1 = slope*ana_stretch(:,1);
hold on;
plot(ana_stretch(:,1),yCalc1,'-k');
hold off;

