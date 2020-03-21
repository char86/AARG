function [NewFigPos] = roiDisplayPosition(screenWidth,screenHeight)

targetPercentWidth = 0.50; targetPercentHeight = 0.6;

targetWidth = targetPercentWidth*screenWidth; targetHeight = targetPercentHeight*screenHeight;

leftPos = (screenWidth*0.5)-(targetWidth*0.5); 
bottomPos = (screenHeight*0.5)-(targetHeight*0.5); 

NewFigPos = [leftPos, bottomPos, targetWidth, targetHeight];