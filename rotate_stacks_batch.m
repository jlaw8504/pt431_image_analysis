function rotate_stacks_batch( directory, centerpoint )
%ROTATE_STACKS_BATCH Creates kymographs from tif stacks using the output
%from population_GUI_v1_2.m
%
%   ROTATE_STACKS_BATCH('directory_name');
%       The function will read in the data needed from the .mat file.
%       Find thec correct image stacks.  Center and rotate the stacks
%       based on the centroid of the plasmid, and the angle of the spindle
%       pole bodies. The name of the .mat files needs to match the base
%       name of the GFP and RFP stacks for batch process to work. Function
%       will search for files only in provided directory.

if nargin == 0
    error('Please give a directory');
elseif nargin == 1
    centerpoint = 'pole';
end

%% Find mat files and corresponding GFP and RFP images
%get files in all folders
files = dir(directory);
%loop through all files
for n = 1:size(files,1)
    %if files are not directories
    if files(n).isdir == 0
        % get the filename
        filename = files(n).name;
        % extension to match
        string = 'mat';
        % split up the string based on .
        split_cell = strsplit(filename,'.');
        % get extension
        extension = split_cell{end};
        % get matching result
        TF_match = strcmp(string,extension);
        % if you are a mat file
        if TF_match == 1
            % load in the data_cell (only variable)
            load(files(n).name,'data_cell');
            % create the names of the TIF stacks that should exist
            string_GFP = strcat(split_cell{1},'_GFP.tif');
            string_RFP = strcat(split_cell{1},'_RFP.tif');
            % loop through files and see if they exist
            for i =1:size(files,1)
                match_name = files(i).name;
                TF_GFP_match = strcmp(string_GFP,match_name);
                if TF_GFP_match == 1
                    GFP_stack = files(i).name;
                end
                TF_RFP_match = strcmp(string_RFP,match_name);
                if TF_RFP_match == 1
                    RFP_stack = files(i).name;
                end
            end
            % check to make sure that both image stacks are in folder
            if exist('GFP_stack','var') == 0
                error(strcat('Could not locate:',string_GFP));
            end
            if exist('RFP_stack','var') == 0
                error(strcat('Could not locate:',string_RFP));
            end
            % read in the TIFF stacks
            GFP_im = tif3Dread(GFP_stack);
            RFP_im = tif3Dread(RFP_stack);
            % change im to double
            GFP_im = im2double(GFP_im);
            RFP_im = im2double(RFP_im);
            % clear for next loop
            clear GFP_stack;
            clear RFP_stack;
            %% Loop through data_cell and alter and save images
            for j=2:size(data_cell,1)
                % Check that both SPBs exist, if not write in blank image
                if isempty(data_cell{j,7}) == 1 ||...
                        isempty(data_cell{j,8}) == 1
                    GFP_im_rot_shift = zeros(size(data_cell{j,12},1),...
                        size(data_cell{j,12},2));
                    RFP_im_rot_shift = zeros(size(data_cell{j,12},1),...
                        size(data_cell{j,12},2));
                else %if two SPB entries exist, do this
                    %% Gather data for rotation
                    % get the correct image plane from stack
                    GFP_im_plane = GFP_im(:,:,data_cell{j,1});
                    % calculate the angle between the SPBs
                    SPB_diff = data_cell{j,7} - data_cell{j,8};
                    angle = atan2d(SPB_diff(2),SPB_diff(1));
                    %% center the images
                    % get image row and col
                    [im_row, im_col] = size(GFP_im_plane);
                    % get centroid from data_cell
                    centroid = data_cell{j,5};
                    % check what centerpoint variable is
                    pole_TF = strcmpi(centerpoint,'pole');
                    plasmid_TF = strcmpi(centerpoint,'plasmid');
                    if plasmid_TF ==1
                        % calc the rows and columns to shift by
                        rows_to_shift = round(im_row/2 - centroid(2));
                        cols_to_shift = round(im_col/2 - centroid(1));
                    elseif pole_TF == 1
                        rows_to_shift = round(im_row/2 - data_cell{j,8}(2));
                        cols_to_shift = round(im_col/2 - data_cell{j,8}(1));
                    end
                    % shift the GFP image
                    GFP_im_shift = circshift(GFP_im_plane,...
                        [rows_to_shift, cols_to_shift]);
                    %% rotate GFP image
                    GFP_im_rot_shift = imrotate(GFP_im_shift,angle,'crop');
                    %% Set up RFP plane
                    % find the RFP plane to use
                    % if the SPBs are in same plane
                    if data_cell{j,7}(3) - data_cell{j,8}(3) == 0
                        % use the plane
                        RFP_im_plane = RFP_im(:,:,data_cell{j,7}(3));
                        % shift the plane
                        RFP_im_shift = circshift(RFP_im_plane,...
                            [rows_to_shift, cols_to_shift]);
                    else
                        % if they are not in same plane, average the image
                        RFP_im_plane_1 = RFP_im(:,:,data_cell{j,7}(3));
                        RFP_im_plane_2 = RFP_im(:,:,data_cell{j,8}(3));
                        RFP_im_plane = (RFP_im_plane_1 + RFP_im_plane_2)/2;
                        % shift the plane
                        RFP_im_shift = circshift(RFP_im_plane,...
                            [rows_to_shift, cols_to_shift]);
                    end % if statement for SPBs in same plane
                    %% Rotate RFP image
                    RFP_im_rot_shift = imrotate(RFP_im_shift,angle,'crop');
                end % if statement for SPBs existing
                %% Add image to a 3D matrix indexed by j-1
                RFP_stack_post(:,:,j-1) = RFP_im_rot_shift;
                GFP_stack_post(:,:,j-1) = GFP_im_rot_shift;
            end
            %% Correct for GFP post stack
            %shift min to 0
            if min(GFP_stack_post(:)) < 0
                GFP_stack_post = GFP_stack_post + abs(min(GFP_stack_post(:)));
            elseif min(GFP_stack_post(:)) > 0 == 1
                GFP_stack_post = GFP_stack_post - min(GFP_stack_post(:));
            end
            %scale max to one
            GFP_stack_post = GFP_stack_post/max(GFP_stack_post(:));
            %convert to unsigned integer 16
            GFP_stack_post = im2uint16(GFP_stack_post);
            
            
            %% Correct for RFP post stack
            %shift min to 0
            if min(RFP_stack_post(:)) < 0
                RFP_stack_post = RFP_stack_post + abs(min(RFP_stack_post(:)));
            elseif min(RFP_stack_post(:)) > 0 == 1
                RFP_stack_post = RFP_stack_post - min(RFP_stack_post(:));
            end
            %scale max to one
            RFP_stack_post = RFP_stack_post/max(RFP_stack_post(:));
            %convert to unsigned integer 16
            RFP_stack_post = im2uint16(RFP_stack_post);
            %% write out the stack
            %Convert post stacks to TIFFs
            for k = 1:size(RFP_stack_post,3)
                image = RFP_stack_post(:,:,k);
                imwrite(image,strcat(split_cell{1},'_RFP_rotated','.tif')...
                    ,'tif','Compression','lzw','WriteMode', 'append');
            end
            clear RFP_stack_post;
            for m = 1:size(GFP_stack_post,3)
                image_2 = GFP_stack_post(:,:,m);
                imwrite(image_2,strcat(split_cell{1},'_GFP_rotated','.tif')...
                    ,'tif','Compression','lzw','WriteMode', 'append');
            end
            clear GFP_stack_post;
        end
    end
end
end


