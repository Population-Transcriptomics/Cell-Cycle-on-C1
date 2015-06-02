macro "SetArea_Measure"{
run("Bio-Formats Macro Extensions"); 
directory=getDirectory("Select a Directory");

//Get the list of all the files in that directory
files = getFileList(directory);

//Start looping through the file names one by one
for (i = 0; i < lengthOf(files); i++) {
		//If the file is bright field image ("BF"), execute the commands in the following loop
		end = endsWith(files[i], "BF.C01");
		if(end == true){
			//Open bright field image file
			//replace "/" with "\\" if you use the macro on windows
			Ext.openImagePlus(directory + "/" + files[i]);
			
			//zoom in
			makeRectangle(125, 125, 250, 250);
			run("To Selection");
				
			//Prepare selection: a circle with default diameter
			defaultDiameter = 13; 
			makeOval(250, 250, defaultDiameter, defaultDiameter);
			
			//Ask the user to select cell location. The size of the circle should not change
			waitForUser("Select cell location with pre-defined circle");
			getSelectionBounds(x, y, wid, hei);
			row = nResults;
			setResult("cell_id", row,substring(files[i], 0, 16));
			setResult("cell_coord_x", row, x);
			setResult("cell_coord_y", row, y);
			close();
			
			//Open the file in green channel for analysis
			//replace "/" with "\\" if you use the macro on windows
			Ext.openImagePlus(directory + "/" + substring(files[i], 0, 16) + "_Green.C01");
			
			//Measure the fluorescence at the location of the cell
			makeRectangle(125, 125, 250, 250);
			run("To Selection");	
			makeOval(x, y, defaultDiameter, defaultDiameter);
			waitForUser("Confirm");
			getStatistics(area, mean, min, max, std);
			setResult("area", row, area);
			setResult("mean_ch2", row, mean);
			setResult("std_ch2", row, std);
			setResult("min_ch2", row, min);
			setResult("max_ch2", row, max);
			
			//Measure the background (bg) fluorescence in a different location
			makeRectangle(125, 125, 250, 250);
			run("To Selection");	
			makeOval(x-50, y+50, defaultDiameter, defaultDiameter);
			waitForUser("Confirm");
			getSelectionBounds(s,t, wid, hei);
			setResult("bg_coord_ch2_x", row, s);
			setResult("bg_coord_ch2_y", row, t);
			getStatistics(area, mean, min, max, std);
			setResult("bg_area_ch2", row, area);
			setResult("bg_mean_ch2", row, mean);
			setResult("bg_std_ch2", row, std);
			setResult("bg_min_ch2", row, min);
			setResult("bg_max_ch2", row, max); 
			close();
			
			//Repeat the measurement of signal and background in the Red channel image
			//replace "/" with "\\" if you use the macro on windows
			Ext.openImagePlus(directory + "/" + substring(files[i], 0, 16) + "_Red.C01");
			makeRectangle(125, 125, 250, 250);
			run("To Selection");
			makeOval(x, y, defaultDiameter, defaultDiameter);
			waitForUser("Confirm");
			getStatistics(area, mean, min, max, std);
			setResult("area_h3", row, area);
			setResult("mean_ch3", row, mean);
			setResult("std_ch3", row, std);
			setResult("min_ch3", row, min);
			setResult("max_ch3", row, max);
			makeRectangle(125, 125, 250, 250);
			run("To Selection");
			makeOval(x-50, y+50, defaultDiameter, defaultDiameter);
			waitForUser("Confirm");
			getSelectionBounds(s, t, wid, hei);
			setResult("bg_coord_ch3_x", row, s);
			setResult("bg_coord_ch3_y", row, t);
			getStatistics(area, mean, min, max, std);
			setResult("bg_area_ch3", row, area);
			setResult("bg_mean_ch3", row, mean);
			setResult("bg_std_ch3", row, std);
			setResult("bg_min_ch3", row, min);
			setResult("bg_max_ch3", row, max); 
			close(); 

			//Ask the user for comment and error report, showing brightfield again
			Ext.openImagePlus(directory + "/" + files[i]);
			title = "No comment";
			width = 512; height = 512;
			Dialog.create("Error report");
			Dialog.addString("Comment:", title);
			Dialog.addChoice("Error:", newArray("0-No Error", "1-No cell", "2-Debris", "3-OutOfFocus","4-MultipleCells"));
  			Dialog.show();
			title = Dialog.getString();
			error = Dialog.getChoice();
			setResult("Error", row, error);
			setResult("Comment", row, title);
			updateResults();
			close();
		}
}
}


