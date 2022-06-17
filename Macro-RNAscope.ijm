//prepare
title=getTitle();
run("Duplicate...", "title=dupl");
run("Set Measurements...", "area limit display redirect=None decimal=3");
run("Set Scale...", "distance=0 known=0 unit=pixel");
roiManager("reset");
run("Split Channels");

//measure DAPI
run("Gaussian Blur...", "sigma=1");
setAutoThreshold("Default dark");
run("Threshold...");
setThreshold(10, 255);
waitForUser("Please check threshold for DAPI");
rename(title+" DAPI area");
run("Measure");
run("Create Selection");
roiManager("Add");
roiManager("Select", 0);
roiManager("Rename", "DAPI");
close();
close();


//measure IL1B
run("Gaussian Blur...", "sigma=1");
setAutoThreshold("Default dark");
run("Threshold...");
setThreshold(20, 255);
waitForUser("Please check threshold for IL1B");
rename(title+" IL1B area");
run("Measure");
run("Create Selection");
roiManager("Add");
roiManager("Select", 1);
roiManager("Rename", "IL1B");
close();

//display
roiManager("Show All");
roiManager("Show None");
roiManager("Show All");
selectWindow("Results");