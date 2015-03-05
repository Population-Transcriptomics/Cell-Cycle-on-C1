macro "SetArea_Measure"{
directory=getDirectory("Select a Directory");
//Get the list of all the files in that direcotry
files=getFileList(directory);
//Start looping through the file names one by one. 
for (i=0; i<lengthOf(files); i++) {
		//If the file is bright field image ("BF"), exectute the commands in the following loop
		end=endsWith(files[i],"BF.C01");
		if(end==true){
			//Open bright field file
			open(directory+"\\"+files[i]);
			//zoom in
			makeRectangle(125, 125, 250, 250);
			run("To Selection");	
			//Prepare selection: a circle with a diameter of 13
			makeOval(250, 250, 13, 13);
			//Ask the user to select cell location. The size of the oval is not changed
			waitForUser("Select cell location with pre-defined circle");
			getSelectionBounds(x,y,wid,hei);
			row=nResults;
			setResult("Label",row,substring(files[i],0,16));
			setResult("Cellcoord.X",row,x);
			setResult("Cellcoord.Y",row,y);
			close();
			//Open the file in green channel for analysis
			open(directory+"\\"+substring(files[i],0,16)+"_Green.C01");
			//Measure the fluorescence at the location of the cell
			makeRectangle(125, 125, 250, 250);
			run("To Selection");	
			makeOval(x, y, 13, 13);
			waitForUser("Confirm");
			getStatistics(area, mean, min, max, std);
			setResult("area", row, area);
			setResult("mean.ch2", row, mean);
			setResult("std.ch2", row, std);
			setResult("min.ch2", row, min);
			setResult("max.ch2", row, max);
			//Measure the background fluorescence in a different location
			makeRectangle(125, 125, 250, 250);
			run("To Selection");	
			makeOval(x-50, y+50, 13, 13);
			waitForUser("Confirm");
			getSelectionBounds(s,t,wid,hei);
			setResult("BGcoord.ch2.X",row,s);
			setResult("BGcoord.ch2.Y",row,t);
			getStatistics(area,mean,min,max,std);
			setResult("bg.area.ch2", row, area);
			setResult("bg.mean.ch2", row, mean);
			setResult("bg.std.ch2", row, std);
			setResult("bg.min.ch2", row, min);
			setResult("bg.max.ch2", row, max); 
			close();
			//Repeat the measurements of signal and background in the Red channel
			open(directory+"\\"+substring(files[i],0,16)+"_Red.C01");
			makeRectangle(125, 125, 250, 250);
			run("To Selection");
			makeOval(x, y, 13, 13);
			waitForUser("Confirm");
			getStatistics(area,mean, min, max, std);
			setResult("area.ch3", row, area);
			setResult("mean.ch3", row, mean);
			setResult("std.ch3", row, std);
			setResult("min.ch3", row, min);
			setResult("max.ch3", row, max);
			makeRectangle(125, 125, 250, 250);
			run("To Selection");
			makeOval(x-50, y+50, 13, 13);
			waitForUser("Confirm");
			getSelectionBounds(s,t,wid,hei);
			setResult("BGcoord.ch3.X",row,s);
			setResult("BGcoord.ch3.Y",row,t);
			getStatistics(area,mean,min,max,std);
			setResult("bg.area.ch3", row, area);
			setResult("bg.mean.ch3", row, mean);
			setResult("bg.std.ch3", row, std);
			setResult("bg.min.ch3", row, min);
			setResult("bg.max.ch3", row, max); 
			//Ask the user for comment and error report
			title = "No comment";
			width=512; height=512;
			Dialog.create("Error report");
			Dialog.addString("Comment:", title);
			Dialog.addChoice("Error:", newArray("0-No Error", "1-No cell", "2-Debris", "3-OutOfFocus","4-MultipleCells"));
  			Dialog.show();
			title = Dialog.getString();
			error = Dialog.getChoice();
			setResult("Error",row,error);
			setResult("Comment", row, title);
			updateResults();
			close();
		}
}
}


