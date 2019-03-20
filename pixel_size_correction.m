%% get the mat files we want to correct
folder = uigetdir('Z:\Josh');
files = dir(folder);
counter = 1;
for n = 1:size(files,1)
    %skip the directories
    if files(n).isdir == 0
        %split up the name to get the file extension
        file_name = files(n).name;
        string_cell = strsplit(file_name,'.');
        extension = string_cell{end};
        %compare string to 'mat'
        TF = strcmp(extension,'mat');
        if TF == 1
            load(files(n).name);
            %% correct the data
            % since there are blanks in the data we need to loop through
            % the cell and test if the cell is actually empty before
            % correcting
            for i = 2:size(data_cell,1)
                %loop through major and minor axis measurements and correct
                data_cell{i,2} = data_cell{i,2} * (64/64.5);
                data_cell{i,3} = data_cell{i,3} * (64/64.5);
                %test if there is a 3D spindle value
                old_spindle = data_cell{i,6};
                if isempty(old_spindle) == 0
                    %loop through and recalcuate the 3D spindle length
                    %get spb1 array
                    spb1 = data_cell{i,7};
                    %convert from pixels to nm using correct pixel size
                    spb1(:,1:2) = spb1(:,1:2)*64;
                    %convert from plane number to nm using 300 nm step size
                    spb1(:,3) = spb1(:,3)*300;
                    %repeat for second spb
                    spb2 = data_cell{i,8};
                    spb2(:,1:2) = spb2(:,1:2)*64;
                    spb2(:,3) = spb2(:,3)*300;
                    %calculate spindle
                    spindle_3D = sqrt((spb1(:,1)-spb2(:,1))^2 + (spb1(:,2)-spb2(:,2))^2 + (spb1(:,3)-spb2(:,3))^2);
                    %replace value
                    data_cell{i,6} = spindle_3D;
                end
                %correct the kMT pixel size
                kmt1 =data_cell{i,9};
                if isempty(kmt1) == 0
                    %square the former answer
                    kmt1_sq = kmt1^2;
                    %find the plane difference from the plasmid and SPB1
                    z = data_cell{i,7}(1,3) - data_cell{i,1};
                    %correct by proper step size (typically 300 nm)
                    z_nm = z * 300;
                    %square the Z component
                    z_nm_sq = z_nm^2;
                    %subtract z_nm_sq from kmt1_sq to get X and Y component
                    kmt1_sq_XY = kmt1_sq - z_nm_sq;
                    %correct the pixel size
                    kmt1_sq_XY_new = ((64^2)/(64.5^2))*kmt1_sq_XY;
                    %add back the Z component
                    kmt1_sq_new = kmt1_sq_XY_new + z_nm_sq;
                    %take sqrt of the X Y and Z component
                    kmt1_new = sqrt(kmt1_sq_new);
                    %assign new kmt1 to data cell
                    data_cell{i,9} = kmt1_new;
                end
                kmt2 = data_cell{i,10};
                if isempty(kmt2) == 0
                    %square the former answer
                    kmt2_sq = kmt2^2;
                    %find the plane difference from the plasmid and SPB1
                    z_2 = data_cell{i,8}(1,3) - data_cell{i,1};
                    %correct by proper step size (typically 300 nm)
                    z_2_nm = z_2 * 300;
                    %square the Z component
                    z_2_nm_sq = z_2_nm^2;
                    %subtract z_nm_sq from kmt1_sq to get X and Y component
                    kmt2_sq_XY = kmt2_sq - z_2_nm_sq;
                    %correct the pixel size
                    kmt2_sq_XY_new = ((64^2)/(64.5^2))*kmt2_sq_XY;
                    %add back the Z component
                    kmt2_sq_new = kmt2_sq_XY_new + z_2_nm_sq;
                    %take sqrt of the X Y and Z component
                    kmt2_new = sqrt(kmt2_sq_new);
                    %assign new kmt1 to data cell
                    data_cell{i,10} = kmt2_new;
                end
            end
            new_name = strcat('corrected_','00',num2str(counter));
            save(new_name,'data_cell');
            counter = counter + 1;
        end
    end
end




