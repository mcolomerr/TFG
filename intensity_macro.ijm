//------------------------------------------------------------------
// FIJI Macro for analysing ROI intensity over time
// First you create a ROI and then the macro will go through the stack
// while you can move it around
// @MColomerRosell
//------------------------------------------------------------------

//Settings of the mesasurements
run("Set Measurements...", "area mean min display redirect=None decimal=3");

//Set area of study
setTool("Rectangle");
waitForUser("SELECT THE AREA OF STUDY");
Roi.getBounds(x, y, width, height);

diameter = 10

//Create a ROI
makeOval(x, y, diameter, diameter);
waitForUser("Create a ROI and move it to the place where the track will start");

interval = 1000

file_name = ""
file_path = "/Volumes/SLN_ICFO_Mariona/ANALYSIS/SMALL_TRACKS/"

//Measure for all the stack
for (i=1; i<(nSlices+1); i++) {
	setSlice(i);
	wait(interval);
	run("Measure");
    }
