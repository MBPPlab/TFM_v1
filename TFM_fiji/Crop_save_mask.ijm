



/* This macro open multicolor  3D files (separate colors) merges the colors and helps croping the cells of interest
It asks for an input and ouput directory
It opens the ROI manager and ask to save in the ROI manager several rectangles corresponding to the cells of interest
Select the regions of interest (ROI) with the “Rectangle” instrument and add them to the ROI list using the key ‘t’. NOTE: when cropping the cell be sure to include a region of at least 5-10 pixel of not moving beads. Exclude  from the analysis cells that are too close to the border or to other cells.
When finish click on ‘OK’.
The macro proposes a mask of the cell: you are satisfied with it click on “OK”. If you are not satisfied click on “Not ok”, select manually a close region with any instrument for selection (e.g. “Freehand” or “Oval”) and click on “Continue”.
*/


// Choose file
path = File.openDialog("Choose a file (multichannel)");
name=File.getName(path);
nameDir=File.getParent(path);
open(path);
Stack.setDisplayMode("grayscale");

name=getTitle();

dirout = getDirectory("Select an output directory ");
name0 = name + '_cell';
nameBeads="movieBeads";
nameMask="movieMask";
print(name0);
// attribute channels correctly
Dialog.create("Attribute channels");
Dialog.addChoice("Beads channel:", newArray("1", "2","3"), "2");
chBeads=parseInt(Dialog.getChoice());
Dialog.addChoice("Mask channel:", newArray("1", "2","3"), "3");
chMask=parseInt(Dialog.getChoice());
//Dialog.addRadioButtonGroup("Type of mask", newArray("Fluo", "Trans"), 1, 2, "Fluo");
Dialog.addCheckbox("Mask channels is fluo (untick if is trans)", 1);
chSel=Dialog.getCheckbox();
Dialog.show();
print(chMask);
print(chSel);

// split channels
selectWindow(name);
run("Duplicate...", "title="+nameBeads+" duplicate channels="+chBeads+" slices=1-numZ");
selectWindow(name);
run("Duplicate...", "title="+nameMask+" duplicate channels="+chMask+" slices=1-numZ");


run("ROI Manager...");
selectWindow(name)
setTool("rectangle");
waitForUser("Select ROIs by clicking on the key T, when ready click OK");



// crop cells
numZ=nSlices;
narea=roiManager("count");
for (i=0; i<narea; i++){
	// all file
	selectWindow(name);
	roiManager("Select", i);
	name1=name0+"_"+i+".tif";
	run("Duplicate...", "title="+name1+" duplicate channels=1-4 slices=1-numZ");
	saveAs("Tiff",dirout+name1);
	selectWindow(name1);
	run("Close");
	// beads
	selectWindow(nameBeads);
	roiManager("Select", i);
	name1=name0+"_"+i+"_movie.tif";
	run("Duplicate...", "title="+name1+" duplicate slices=1-numZ");
	run("Channels Tool...");
	saveAs("Tiff",dirout+name1);
	selectWindow(name1);
	run("Close");
	// mask
	selectWindow(nameMask);
	roiManager("Select", i);
	name1=name0+"_"+i+"_Mask.tif";
	run("Duplicate...", "title="+name1+" duplicate slices=1-numZ");
	run("Channels Tool...");
	saveAs("Tiff",dirout+name1);
	selectWindow(name1);
	run("Close");
}
selectWindow(nameMask);
run("Close");
selectWindow(nameBeads);
run("Close");


ch=getBoolean("Which file are you using for the mask?","Trans", "Fluo");


// create mask cell by cell
for (i=0; i<narea; i++){
	name1=name0+"_"+i+"_Mask.tif";
	open(name1);
	selectWindow(name1);
//	roiManager("Select", i);
//	run("Duplicate...", "title=tmp"+file_for_mask+" duplicate");	
	if (chSel==0){
		mask_trans(name1);}
	else {
		mask_fluo(name1);
	}
	a=getBoolean("Are you happy with the mask?");
	if (a==0){
		close(name1);
		open(name1);
		waitForUser("Select the ROI manually and click OK (it will be the same for all timepoints)");
		setTool("oval");
		setBackgroundColor(0, 0, 0);
		setForegroundColor(255, 255, 255);
		run("Clear Outside", "stack");
		run("Fill", "stack");
		run("Make Binary", "method=Default background=Default calculate black");
	}
	saveAs("Tiff", dirout+name1);
	run("Close");	
}


function conversion8bit(){
	run("Enhance Contrast...", "saturated=0 normalize process_all use");
	run("Divide...", "value=256 stack");
	run("Enhance Contrast...", "saturated=0 normalize process_all use");
	run("8-bit");

}	

function mask_trans(file_for_mask){	
getDimensions(width, height, channels, slices, frames);
if (frames<slices){
	Stack.setDimensions(channels, frames, slices);
	getDimensions(width, height, channels, slices, frames);
}

// work on trans
run("Variance...", "radius=5 stack");
run("Gaussian Blur...", "sigma=3 stack");
run("Auto Threshold", "method=Otsu white stack");
run("Dilate", "stack");
run("Close-", "stack");
run("Fill Holes", "stack");
run("Erode","stack");
run("Analyze Particles...", "size=1000-Infinity pixel clear add stack");
for (i=1; i<frames+1; i++){
	roiManager("Select", i-1);
	//Stack.setFrame(i); 
	run("Clear Outside", "slice");	
}
roiManager("Show All");
roiManager("Show None");
selectWindow("ROI Manager");
run("Close");
selectWindow(file_for_mask);
}



function mask_fluo(file_for_mask){	
getDimensions(width, height, channels, slices, frames);
if (frames<slices){
	Stack.setDimensions(channels, frames, slices);
	getDimensions(width, height, channels, slices, frames);
}

// work on fluo channel
run("Gaussian Blur...", "sigma=3 stack");
run("Auto Threshold", "method=Otsu white stack");
run("Dilate", "stack");
run("Close-", "stack");
run("Fill Holes", "stack");
run("Erode","stack");
run("Analyze Particles...", "size=1000-Infinity pixel clear add stack");
for (i=1; i<frames+1; i++){
	roiManager("Select", i-1);
	//Stack.setFrame(i); 
	run("Clear Outside", "slice");	
}
roiManager("Show All");
roiManager("Show None");
selectWindow("ROI Manager");
run("Close");
selectWindow(file_for_mask);
}







