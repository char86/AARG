function [GUI] = ResizeGUIForCurrentScreen(GUI, screenWidth, screenHeight, widthScaleExtension, heightScaleExtension)

% Outline
% ResizeGUIForCurrentScreen attempts to fit the GUI to the current screen.
% AARG development work was done mostly on a screen with screenWidth = 1440
% pixels and screenHeight = 900 pixels. This gives a ratio of width/height
% ratio of 1.6. If the ratio is greater than this widthScaleExtension = 1.3
% and heightScaleExtension = 1.1. A catelogue of different combinations of
% extension values might need to be developed for different ratios. 


standardHeight = 900; standardWidth = 1440;

whRatio = screenWidth/screenHeight; approxStandardWidth = standardHeight*whRatio;
factorW = approxStandardWidth/standardWidth;

widthScale = ((screenWidth/approxStandardWidth)*factorW)*widthScaleExtension; 
heightScale = (screenHeight/standardHeight)*heightScaleExtension;

guiChildren = GUI.Children; numberOfGUIChildren = size(GUI.Children,1);

oldGUIPosition = GUI.Position;
GUI.Position = [oldGUIPosition(1)*widthScale oldGUIPosition(2)*heightScale oldGUIPosition(3)*widthScale oldGUIPosition(4)*heightScale];

for childIdxF1 = 1 : numberOfGUIChildren
    oldChildPosition = guiChildren(childIdxF1).Position;
    guiChildren(childIdxF1).Position = [oldChildPosition(1)*widthScale,... 
                                      oldChildPosition(2)*heightScale,...
                                      oldChildPosition(3)*widthScale,...
                                      oldChildPosition(4)*heightScale];
    f1Children = guiChildren(childIdxF1).Children;
    numberOfF1Children = size(guiChildren(childIdxF1).Children, 1);
    
    try 
        if factorW < 1
            guiChildren(childIdxF1).FontSize = 10;
        end
    catch
    end

    for childIdxF2 = 1 : numberOfF1Children
        oldF1ChildPosition = f1Children(childIdxF2).Position;
        f1Children(childIdxF2).Position = [oldF1ChildPosition(1)*widthScale,...
                                           oldF1ChildPosition(2)*heightScale,...
                                           oldF1ChildPosition(3)*widthScale,...
                                           oldF1ChildPosition(4)*heightScale];
                                       
        try 
            if factorW < 1
                f1Children(childIdxF2).FontSize = 10;
            end
        catch
        end                          
                                       
    end
    guiChildren(childIdxF1).Children = f1Children;
end
GUI.Children = guiChildren;
