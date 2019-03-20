function rotate_center_image(input_cell)
%%This function will take the data_cell cell array of the population_gui
%%and output a TIF stack of centered and rotated images of the plasmid

%% Gather parameters for the FOR LOOP
%get the row size of the cell array
cell_row = size(input_cell,1);
%get the size fo the image
[im_row, im_col] = size(input_cell{2,12});
%pre-allocate the size of the 3D matrix to which the images will go
im_mat = zeros(im_row,im_col,(cell_row -1));
%% For loop starts here
% since row 1 contains labels start at n = 2
for n = 2:cell_row
    %% shift cell to center
    % read in the image
    im = input_cell{n,12};
    %get centroid from input_cell
    centroid = input_cell{n,5};
    %calc the rows and columns to shift by
    rows_to_shift = round(im_row/2 - centroid(2));
    cols_to_shift = round(im_col/2 - centroid(1));
    im_shift = circshift(im, [rows_to_shift, cols_to_shift]);
    %% calc the angle of the plasmid to the X-axis
    %create a binary image
    im_bin = im_shift > 0;
    %check that there is only one object
    CC = bwconncomp(im_bin);
    num_obj = CC.NumObjects;
    %create the morphological operator to erode and dilate with
    SE_disk = strel('disk', 1);
    if num_obj > 1
        im_bin = imerode(im_bin, SE_disk);
        im_bin = imdilate(im_bin, SE_disk);
        im_shift = im_shift .* im_bin;
    else
        im_bin = imdilate(im_bin, SE_disk);
        im_bin = imerode(im_bin, SE_disk);
        im_shift = im_shift .* im_bin;
    end
    %use regionprops to get angle
    S = regionprops(im_bin,'Orientation');
    angle = S.Orientation;
    im_rot = imrotate(im_shift,(-angle),'crop');
    %push the im_rot to the proper position in the matrix
    im_mat(:,:,(n-1)) = im_rot;
end
%% Convert matrix into a TIFF stack
%add in a backround to make image similar to simulated image output
% im_mat = im_mat + randi(round(max(im_mat(:))),im_row,im_col,(cell_row-1));
im_mat = im_mat + max(im_mat(:))/2;
%set first voxel to 0 so shift min (next code block) doesn't git rid of the
%effect
im_mat(1,1,:) = 0;

%shift min to 0
if min(im_mat(:)) < 0
    im_mat = im_mat + abs(min(im_mat(:)));
elseif min(im_mat(:)) > 0 == 1
    im_mat = im_mat - min(im_mat(:));
end
%scale max to one
im_mat = im_mat/max(im_mat(:));

%convert to unsigned integer 16
im_mat = im2uint16(im_mat);

%Convert to tif and save
[filename, pathname] = uiputfile('*.tif', 'Where would you like to save the file?');
for i = 1:size(im_mat,3)
    image = im_mat(:,:,i);
    imwrite(image,strcat(pathname,filename),'tif','Compression',...
        'lzw','WriteMode', 'append');
end
end


