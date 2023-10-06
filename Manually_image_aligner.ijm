// Prompt the user to select a directory
inputDir = getDirectory("Choose a Directory");

// Number of channels
n_channels = 4;

// Get a list of image file names in the selected directory
list = getFileList(inputDir);

// Sort the list alphabetically
Array.sort(list);

// Initialize variables to store the biggest dimensions
biggest_x = 0;
biggest_y = 0;

// Loop through each image file again to perform operations
for (i = 0; i < list.length; i++) {
    if (endsWith(list[i], ".tif") || endsWith(list[i], ".jpg") || endsWith(list[i], ".png")) {
        filePath = inputDir + "/" + list[i];
        open(filePath);
        
        // Set the tool to "multipoint"
        setTool("multipoint");
        
        // Wait for the user to make points
        waitForUser("Make points on the image. Press OK when done.");
        
        // Get the list of points
        getSelectionCoordinates(xCoordinates, yCoordinates);
        
        if (xCoordinates.length > 0) {
            // Calculate translation required to center the selected point
            centerX = Math.round(getWidth() / 2);
            centerY = Math.round(getHeight() / 2);
            translateX = Math.round(centerX) - xCoordinates[0];
            translateY = Math.round(centerY) - yCoordinates[0];
                        
            // Calculate the new canvas size
            newWidth = getWidth() + Math.abs(translateX) * 2;
            newHeight = getHeight() + Math.abs(translateY) * 2;
           
            // Resize the canvas
            run("Canvas Size...", "width=" + newWidth + " height=" + newHeight + " position=Center zero");
            
            // Translate the image
            run("Translate...", "x=" + translateX + " y=" + translateY + " interpolation=None");
            

            
            // Compare with biggest dimensions and update if needed
	        if (getWidth() > biggest_x) {
	            biggest_x = getWidth();
	        }
	        if (getHeight() > biggest_y) {
	            biggest_y = getHeight();
	        }
            
            // Overwrite the image
            saveAs("Tiff", filePath);
            close();
        }
        
        // Clear the points
        run("Select None");
    }
}


// Loop through each image file and open it to make the dimentions the same
for (i = 0; i < list.length; i++) {
    if (endsWith(list[i], ".tif") || endsWith(list[i], ".jpg") || endsWith(list[i], ".png")) {
        filePath = inputDir + "/" + list[i];
        open(filePath);
        
        // Resize the canvas
        run("Canvas Size...", "width=" + biggest_x + " height=" + biggest_y + " position=Center zero");
        
        // Overwrite the image
        saveAs("Tiff", filePath);
        close();
    }
}

// Open the sequence
File.openSequence(inputDir, "virtual");
run("Stack to Hyperstack...", "order=xyczt(default) channels=" + n_channels + " slices=" + list.length + " frames=1 display=Color");
Stack.setChannel(1);
run("Blue");
run("Enhance Contrast", "saturated=0.35");
Stack.setChannel(2);
run("Green");
resetMinAndMax();
Stack.setChannel(3);
run("Red");
resetMinAndMax();
Stack.setChannel(4);
run("Magenta");
resetMinAndMax();

// Save final stack
saveAs("Tiff", inputDir + "/stack.tif");
