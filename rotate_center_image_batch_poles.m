function rotate_center_image_batch_poles(directory)
%%This function will take the data_cell cell array of the population_gui
%%and output a TIF stack of centered and rotated images of the plasmid and
%%a TIFF of a maximum intensity, no BG subtraction, kymograph.

%% Loop through directory contents and find the data files
%change to that directory;
cd(directory);
%get the files from that directory
files = dir(directory);
%looking for mat files, so set the string to 'mat'
string = 'mat';
for j=1:size(files,1)
    if files(j).isdir == 0
        string_cell = strsplit(files(j).name, '.');
        string_ext = string_cell(end);
        match_TF = strcmp(string,string_ext);
        if match_TF == 1
            load(files(j).name);
            %% Gather parameters for the FOR LOOP
            % Limit the output to only timelapses that do no split
            % read in the split column
            split_array = cell2mat(data_cell(2:end,11));
            split_check = sum(split_array);
            if split_check == 0
                %get the row size of the cell array
                cell_row = size(data_cell,1);
                %get the size fo the image
                [im_row, im_col] = size(data_cell{2,12});
                %pre-allocate the size of the 3D matrix to which the images will go
                im_mat = zeros(im_row,im_col,(cell_row -1));
                im_mat_poles = im_mat;
                %% For loop starts here
                % since row 1 contains labels start at n = 2
                for n = 2:cell_row
                    %% shift cell to center
                    % read in the image
                    im = data_cell{n,12};
                    %create a SPB only image
                    if isempty(data_cell{n,7}) == 0
                        im_poles = zeros(im_row,im_col);
                        spb1 = data_cell{n,7};
                        im_poles(spb1(2)-1:spb1(2)+1,spb1(1)-1:spb1(1)+1) = max(im(:));
                    end
                    if isempty(data_cell{n,8}) == 0
                        if exist('im_poles','var') == 0
                            im_poles = zeros(im_row,im_col);
                        end
                        spb2 = data_cell{n,8};
                        im_poles(spb2(2)-1:spb2(2)+1,spb2(1)-1:spb2(1)+1) = max(im(:));
                    end
                    if isempty(data_cell{n,7}) == 1 && isempty(data_cell{n,8}) == 1
                        im_poles = zeros(im_row,im_col);
                    end
                    %im is the plasmid only image
                    im = data_cell{n,12};
                    %get centroid from data_cell
                    centroid = data_cell{n,5};
                    %calc the rows and columns to shift by
                    rows_to_shift = round(im_row/2 - centroid(2));
                    cols_to_shift = round(im_col/2 - centroid(1));
                    im_shift_poles = circshift(im_poles, [rows_to_shift, cols_to_shift]);
                    im_shift = circshift(im, [rows_to_shift, cols_to_shift]);
                    %% calc the angle of the plasmid to the X-axis
                    %create a binary image of plasmid only image
                    im_bin = im_shift > 0;
                    %check that there is only one object
                    CC = bwconncomp(im_bin);
                    num_obj = CC.NumObjects;
                    %create the morphological operator to erode and dilate
                    %with structural element on image with poles
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
                    im_rot_poles = imrotate(im_shift_poles,(-angle),'crop');
                    im_rot = imrotate(im_shift,(-angle),'crop');
                    %push the im_rot to the proper position in the matrix
                    im_mat(:,:,(n-1)) = im_rot;
                    im_mat_poles(:,:,(n-1) )= im_rot_poles;
                    %clear im_poles for an above if statement
                    clear im_poles
                end
                %% Convert matrix into a TIFF stack
                %add in a backround to make image similar to simulated image output
                % im_mat = im_mat + randi(round(max(im_mat(:))),im_row,im_col,(cell_row-1));
                %                 im_mat = im_mat + max(im_mat(:))/2;
                %set first voxel to 0 so shift min (next code block) doesn't git rid of the
                %effect
                %                 im_mat(1,1,:) = 0;
                %% Correct for plasmid image
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

                
                %% Correct for pole image
                %shift min to 0
                if min(im_mat_poles(:)) < 0
                    im_mat_poles = im_mat_poles + abs(min(im_mat_poles(:)));
                elseif min(im_mat_poles(:)) > 0 == 1
                    im_mat_poles = im_mat_poles - min(im_mat_poles(:));
                end
                %scale max to one
                im_mat_poles = im_mat_poles/max(im_mat_poles(:));
                %convert to unsigned integer 16
                im_mat_poles = im2uint16(im_mat_poles);
                
                %                 %create the kymograph matrix
                %                 %preallocate kymo_mat for speed
                %                 kymo_mat = zeros(size(im_mat,3),size(im_mat,2));
                %                 for k = 1:size(im_mat,3)
                %                     kymo_mat(k,:) = max(squeeze(im_mat(:,:,k)));
                %                 end
                
                
                
                %Convert im_mat to tif and save
                for i = 1:size(im_mat,3)
                    image = im_mat(:,:,i);
                    imwrite(image,strcat(string_cell{1},'.tif'),'tif','Compression',...
                        'lzw','WriteMode', 'append');
                end
                for i = 1:size(im_mat_poles,3)
                    image_poles = im_mat_poles(:,:,i);
                    imwrite(image_poles,strcat(string_cell{1},'_poles','.tif'),'tif','Compression',...
                        'lzw','WriteMode', 'append');
                end
            end
        end
    end
end

