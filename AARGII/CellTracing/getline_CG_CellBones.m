function [xs,ys, linehandle, linehandle_error, seltype, abortswitch, lastAction] = getline_CG_CellBones(axishandle,dospline,...
    CA_Bones, seltype, StartPointColor,linehandle_error, abortswitch, lastAction, fh)

%============================================================================
% [xs,ys, linehandle] = getpoints(axishandle,dospline)
% select an area interactively with mouse clicking a polygon
% right click to close the selection
% Jianhua Cang, Dec-12-2005, with original code from Matt Caywood
%============================================================================
%CG: This script has been modified by Charlie J Gilbride for use with the
%main function CellBones. Original comments that do not conflict with the
%modifications have been left in place. 
%============================================================================
% Find parent figure for the argument axishandle
%============================================================================
axes(axishandle); figure(get(axishandle, 'Parent'));

hold on
pointhandles = []; xpts = []; ypts = []; outlinehandle = []; n = 0;
but = 1; BUTN = 0; KEYB = 1; done = 0; 
KeyPressRes = '';

%===========================================================================
% Loop until right hand mouse button or keayboard is pressed
%===========================================================================
while ~done && ~abortswitch 
 %===========================================================================
 % Analyze each buttonpressed event
 %==========================================================================
%  try
     if strcmp(lastAction, 'a-press-but-no-alt-click')
         keyb_or_butn = 2; %CG: neither keyboard press, nor mouse click
     elseif ~strcmp(lastAction, 'a-press-but-no-alt-click')
        keyb_or_butn = waitforbuttonpress;
     end
    if ~strcmp(linehandle_error, 'noerror')
%CG: replace title instructons if the last action triggered an error.
        title(sprintf(strcat('Instructions: Left-click on start point or white line to begin tracing a dendrite.\n',...
        'Build line with further left-clicks. Right-click to stop building current white line.\n',...
        'Left-click on start point or white line again to build new white line.\n',...
        'Alternatively, press "a" to accept and move to next step')), 'FontSize', 16)
    end
    PointerStatus = get(gcf, 'Pointer');
    if strcmp(PointerStatus, 'arrow')
        seltype = get(gcf,'SelectionType');
        if keyb_or_butn == BUTN && ~strcmp(lastAction,'a-press-but-no-alt-click');
          switch seltype
          case 'normal',
            but = 1; currpt = get(axishandle, 'CurrentPoint');
            linehandle_error = 'noerror'; lastAction = 'normal';
          case 'alt',
            linehandle_error = 'noerror'; but = 2; done = 1; lastAction = 'alt';
          otherwise,
            but = 2; lastAction = 'other';
          end;           
        elseif keyb_or_butn == KEYB && ~strcmp(lastAction,'a-press-but-no-alt-click')
            KeyPressRes = get(gcf,'CurrentCharacter');
            if strcmp(KeyPressRes, 'a') %&& strcmp(lastAction, 'alt')
%CG: if the last mouse click was right and the last action was an "a" key press
%but=2 
                if strcmp(lastAction, 'alt') || isempty(lastAction)
%CG: if strcmp(lastAction, 'alt') then the user has behaved as expected and closed the 
%last line before pressing 'a'. if isempty(lastAction), this could mean
%that the user has chosen not to add any additional white line. 
                    but = 3;
                elseif strcmp(lastAction, 'normal')
                    but = 1; currpt = get(axishandle, 'CurrentPoint');
                    linehandle_error = 'noerror'; lastAction = 'a-press-but-no-alt-click';
                end
               
            elseif ~isempty(KeyPressRes)
                linehandle_error = 'error2'; linehandle = []; but = 2; done = 1;
            end
        elseif strcmp(lastAction,'a-press-but-no-alt-click')
            but = 3;
        end
         %===========================================================================
         % Get coordinates of the last buttonpressed event
         %===========================================================================
        if but == 1
            xi = round(currpt(2,1)); yi = round(currpt(2,2));
        end
        cObj = gco(fh);
        if (strcmp(get(cObj,'type'), 'line') || strcmp(get(cObj,'type'), 'image')) && but == 1
            Sz_CA_Bones = size(CA_Bones); NumBones = Sz_CA_Bones(1); cBoneCell = 0;
            while cBoneCell < NumBones
                cBoneCell = cBoneCell + 1;
                currentLine_XData = get(CA_Bones{cBoneCell, 1}, 'XData');
                currentLine_YData = get(CA_Bones{cBoneCell, 1}, 'YData');
                
                cOdj_XData = get(cObj, 'XData'); cObj_YData = get(cObj, 'YData');

                if sum(currentLine_XData) - sum(cOdj_XData) == 0 && sum(currentLine_YData) - sum(cObj_YData) == 0
%CG: Make sure the current line is the clicked object
                    
                    AllLineSubs = CA_Bones{cBoneCell, 3}; 
                    Sz_AllLineSubs = size(AllLineSubs); NumPnts = Sz_AllLineSubs(1);

                    currpt_exp = ones(NumPnts,1); currpt_X_exp = currpt_exp.*xi;
                    currpt_Y_exp = currpt_exp.*yi;

                    diff_cp_X = currpt_X_exp - AllLineSubs(:,1);
                    diff_cp_Y = currpt_Y_exp - AllLineSubs(:,2);

                    [~, cloIdx_X] = min(abs(diff_cp_X));
                    [~, cloIdx_Y] = min(abs(diff_cp_Y));
                    currpt
                    if cloIdx_X ~= cloIdx_Y
                        p = randperm(2,1);
%CG: not clear why, but even if you take the coordinate pair with the shortest
%diagonal distance between it and the selected point getline does not
%always start a new white line anywhere near the selected spot on a
%pre-existing white line. 
                        if p == 1
                            xi = AllLineSubs(cloIdx_X,1); yi = AllLineSubs(cloIdx_Y,2);
                        else
                            xi = AllLineSubs(cloIdx_X,1); yi = AllLineSubs(cloIdx_Y,2);
                        end
                    else
%CG: update the current point coordinates so they connect the new line with
%the old one. 
                        xi = AllLineSubs(cloIdx_X,1); yi = AllLineSubs(cloIdx_Y,2);
                        cBoneCell = NumBones;
                    end
                end
            end
        elseif strcmp(get(cObj,'type'), 'patch') && but == 1
            
            cObjColor = get(cObj, 'FaceColor'); ColorTest = cObjColor - StartPointColor;

            if sum(ColorTest) == 0 && isempty(xpts)
                PatchXData = get(cObj, 'XData'); PatchYData = get(cObj, 'YData');

                if numel(PatchXData) > 1
                    xi = PatchXData(4); yi = PatchYData(2);
                elseif numel(PatchXData) == 1
                    xi = PatchXData; yi = PatchYData;
                end
                
            elseif sum(ColorTest) ~= 0 && isempty(xpts)
               linehandle_error = 'error1'; linehandle = [];but = 2; done = 1;
            end

        end
        get(cObj,'type')
        if ~strcmp(get(cObj,'type'), 'patch') &&...
            ~strcmp(get(cObj,'type'), 'line') && isempty(xpts)
            if ~isempty(KeyPressRes)
                if ~strcmp(KeyPressRes,'a')
                    linehandle_error = 'error1'; linehandle = []; but = 2; done = 1;
                end
            elseif isempty(KeyPressRes)
                linehandle_error = 'error1'; linehandle = []; but = 2; done = 1;
            end
        end
         %===========================================================================
         % Start a spline throught the points or
         % update the line through the points with a new spline
         %===========================================================================
         if but ==1
           if ~isempty(outlinehandle)
              delete(outlinehandle);
           end
           pointhandles(n+1) = plot(xi,yi,'wo','MarkerSize',4,'LineWidth',4);
           n = n+1; xpts(n,1) = round(xi); ypts(n,1) = round(yi);
           %===========================================================================
           % Draw a spline line through the points
           %===========================================================================
           if n > 1
             t = 1:n; ts = 1: 0.1 : n;
             if (dospline)
                 xs = spline(t, xpts, ts); ys = spline(t, ypts, ts);
                 outlinehandle = plot(xs,ys,'r-');
             else
                 outlinehandle = line(xpts,ypts,'Color',[1 1 1],'LineWidth',2);
             end
           end;
         elseif but == 3
             %===========================================================================
                 % Exit for "a" keyboard input 
             %===========================================================================
             seltype = 'finish'; linehandle_error = 'noerror'; done = 1;
         end;
    end
%  catch
%      abortswitch = 1;
%  end
end;

if abortswitch == 0
%===========================================================================
% (re)draw the final spline (CG: if error-free)
%===========================================================================
    if strcmp(linehandle_error, 'noerror')
       if ~isempty(outlinehandle); delete(outlinehandle); end;
       linehandle = line(xpts,ypts,'Color',[1 1 1],'LineWidth',2);
       xs = xpts; ys = ypts;
    elseif strcmp(linehandle_error, 'error1') || strcmp(linehandle_error, 'error2')
        xs = []; ys = [];
    end
    drawnow;
%===========================================================================
% Delete the point markers
%===========================================================================
    if ~isempty(pointhandles); delete(pointhandles); end;
elseif abortswitch == 1
    xs = [];ys = []; linehandle = []; linehandle_error = ''; seltype = '';
end








