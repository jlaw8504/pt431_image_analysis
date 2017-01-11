function im_props = pt431_image_analysis(directory)

%This function will parse a directory of MicroscopeSimulator 2 generated
%images, pick the in-focus plane of each stack, and measure the major and
%minor axis of a thresholded image.

%This function uses tif3Dread.m provided by MBL's CIAN course.

%% Parse the directory
% get a list of tif files in directory
tif_files = dir(fullfile(directory,'*.tif'));
for n = 1:length(tif_files)
    %Parse the 3D stack to a 3D matrix of uint16
    im = tif3Dread(fullfile(directory,tif_files(n).name));
    im_dbl = im2double(im);
    %find the third dimension with the brightest pixel
    [~,max_idx] = max(im_dbl(:));
    [~,~,inf_plane] = ind2sub(size(im_dbl),max_idx);
    im_inf = im_dbl(:,:,inf_plane);
    %subtract the mode from the image as BG subtraction
    im_mode = mode(im_inf(:));
    im_bg = im_inf - im_mode;
    %convert any negative numbers to zeros
    im_bg = im_bg .* (im_bg > 0);
    %Determine the plasmid signal by using Otsu threshold
    thresh = multithresh(im_bg);
    im_bin = im_bg > thresh;
    %use regionprops to determine the major/minor axis, centroid and
    %orientation
    im_props(n)= regionprops(im_bin,'MinorAxisLength','MajorAxisLength',...
        'Centroid','Orientation');
    AspectRatio{n} = ...
        im_props(n).MajorAxisLength/im_props(n).MinorAxisLength;
end
    [im_props.AspectRatio] = AspectRatio{:};
end