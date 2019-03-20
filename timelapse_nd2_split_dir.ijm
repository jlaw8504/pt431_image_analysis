dir1 = getDirectory("Choose Source Directory ");
list = getFileList(dir1);
for (n=0; n<list.length; n+=2) {
     run("Bio-Formats Importer", "open="+dir1+list[n]+" autoscale color_mode=Grayscale concatenate_series open_all_series rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
     name = File.nameWithoutExtension;
     nameext = File.name;
     print(n);
     roiManager("Open",dir1+list[n+1]);

//Get the name of the open image
filename = getInfo("image.filename");
//remove the file extension
name_array = split(filename,".");
basename = name_array[0];
//Get the directory of the open image
file_dir = getInfo("image.directory");
//Get the number of objects in the ROI Manager
count = roiManager("count");
//Loop through the ROIs in the ROI Manager
for (i=0; i<count; i++) {
	j = i+1;
	roiManager("Select", i);
	run("Duplicate...", "title=Cell duplicate");
	run("Stack Slicer", "split_channels stack_order=XYCTZ");
	if (i<10) {
		string = "Cell_00";
	} else {
		string = "Cell_0";
	}
	selectWindow("Cell - C=2");
	saveAs("Tiff", file_dir+string+j+"_"+basename+"_RFP");
	close(string+j+"_"+basename+"_RFP"+".tif");
	selectWindow("Cell - C=1");
	saveAs("Tiff", file_dir+string+j+"_"+basename+"_GFP");
	close(string+j+"_"+basename+"_GFP"+".tif");
	selectWindow("Cell - C=0");
	saveAs("Tiff", file_dir+string+j+"_"+basename+"_Trans");
	close(string+j+"_"+basename+"_Trans"+".tif");
}
// delete roi's before loading in more
roiManager("reset");
title = getTitle();
close(title);
}
