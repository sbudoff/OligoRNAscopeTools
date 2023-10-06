// Prompt the user to select a directory
inputDir = getDirectory("Choose a Directory");

// Get a list of image file names in the selected directory
list = getFileList(inputDir);


// Loop through each image file again to perform operations
for (i = 0; i < list.length; i++) {
    if (endsWith(list[i], ".tif") || endsWith(list[i], ".jpg") || endsWith(list[i], ".png")) {
    	in_name = list[i];
    	small_name = in_name.substring(0,in_name.length-4)+"_small.tif";
    	red_name = in_name.substring(0,in_name.length-4)+"_small_rbms.tif";
    	png_name = in_name.substring(0,in_name.length-4)+"_small_rbms.png";
        filePath = inputDir + "/" + in_name;
        open(filePath);
        W = Math.round(getWidth()/10);
        H = Math.round(getHeight()/10);
        run("Scale...", "x=- y=- z=1.0 width=" + W + " height=" + H + " depth=4 interpolation=Bilinear average create");
		selectWindow(in_name);
		close();
		saveAs("Tiff", inputDir + "/" + small_name);
		run("Split Channels");
		selectWindow("C4-"+small_name);
		close();
		selectWindow("C1-"+small_name);
		close();
		selectWindow("C2-"+small_name);
		close();
		selectWindow("C3-"+small_name);
		saveAs("Tiff", inputDir + "/" + red_name);
		run("Enhance Contrast", "saturated=0.35");
		run("RGB Color");
		saveAs("PNG", inputDir + "/" + png_name);
		close();
        
        
    }
}