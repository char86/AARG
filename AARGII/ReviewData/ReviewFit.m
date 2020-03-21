function ReviewFit(cDir,suffixStrs,text_status,cConditionIdx,frameRate)

%Outline: ReviewFit uses a simple gui enabling the user to rapidly check the
%accuracy of curve fitting and peak detection. However, the user also has
%the option to reverse the previous action. 

textleft = 2.1;
text_top = 1.25;
textwidth = 0.3;
textheight = 0.1;
textvertspace = 0.3;

fontsz = 16;

if nargin == 0
    suffixStrs = {'a_baseline', 'b_TTAP2'};
end

%Set up gui:
scrsz=get(0,'ScreenSize'); figpos = [scrsz(3)*0.2,scrsz(4)*0.2,scrsz(3)*0.6,scrsz(4)*0.6];
fh = figure('NumberTitle','off','Visible','off','Position',figpos);

winwidth = figpos(3); winheight = figpos(4);

instructions_textloc = round([(figpos(1)*textleft),...
    figpos(2)*(text_top+0.5),winwidth*(textwidth),winheight*(textheight+0.25)]);

instructions_text = uicontrol('Style','text',...
         'Position',instructions_textloc, 'FontSize', 10,...
         'HorizontalAlignment', 'left');   
str = sprintf(strcat('INSTRUCTIONS:\nPress "a" key to accept trace.\n\n',...
    'Press "v" key to reject only curve fit.\n\n',...
    'Press "r" key to reject all measurements for current trace.\n\n',...
    'Press "shift+u" keys to undo. Press "shift+r" to restart.\n\n',...
    'Pressing "a", "v" or "r" will also prompt presentation of next trace.'));
instructions_text.String = str;


accept_textloc = round([(figpos(1)*textleft),...
    figpos(2)*text_top,winwidth*textwidth,winheight*textheight]);
accept_text = uicontrol('Style','text',...
         'Position',accept_textloc, 'FontSize', fontsz,...
         'HorizontalAlignment', 'left');     

newtext_top = text_top-textvertspace;

acceptamponly_textloc = round([(figpos(1)*textleft),...
    figpos(2)*newtext_top,winwidth*textwidth,winheight*textheight]);
acceptamponly_text = uicontrol('Style','text',...
         'Position',acceptamponly_textloc, 'FontSize', fontsz,...
         'HorizontalAlignment', 'left'); 
    
newtext_top = text_top-(textvertspace*2);

reject_textloc = round([(figpos(1)*textleft),...
    figpos(2)*newtext_top,winwidth*textwidth,winheight*textheight]);
reject_text = uicontrol('Style','text',...
         'Position',reject_textloc, 'FontSize', fontsz,...
         'HorizontalAlignment', 'left'); 
    
axes_loc = round([winwidth*0.075, winheight*0.15, winwidth*0.6, winheight*0.75]);   
ha = axes('Units','pixels','Position',axes_loc);
align([accept_text,acceptamponly_text,reject_text],'Left','None');
%Initiate gui: 
%***Integrate later...
if nargin > 0
    Dirs = {}; Dirs{1} = cDir;
else
    Dirs = uipickfiles('prompt','Review traces');
    % Dirs = {}; Dirs{1} ='/Users/cjgilbride/Documents/DecJan2016_17/MATLAB/AARG_pub_datatest/CG0706173';
end
if iscell(Dirs)
    NumDirs = size(Dirs,2); NumConds = size(suffixStrs,2);
    for cDirIdx = 1 : NumDirs
        
        datadir = Dirs{cDirIdx}; cd(Dirs{cDirIdx}); [~,CellName,~] = fileparts(datadir);
        %cConditionIdx = 0;
        while cConditionIdx <= NumConds 
            
            fh.Name = strcat('Experiment: "',CellName, '"; Condition: "',suffixStrs{cConditionIdx},'"');
            if nargin > 0
                text_status.String = strcat('Reviewing data for experiment: "',...
                    CellName, '"; Condition: "',suffixStrs{cConditionIdx},'"');
            end
            filestr = strcat(CellName,suffixStrs{cConditionIdx},'_MeasuresRAW.mat');
            data = load(filestr); CA_Measure = data.CA_Measure;

            Sz_CA_Measure = size(CA_Measure); totaltraces = 0; 
            roi_indices = []; CA_Segs = {};
            roi_idx = []; trace_idx = 1;
            reviewedtraces = 1; acceptedtraces = 0; partacceptedtraces = 0; rejectedtraces = 0;
        %CG: get number of trace segments the user must inspect.
            try
                
%                 BreakTry
               undostring = data.undostring; parameters = data.parameters;
            catch
                undostring = '';
                %CG: undostring carries a record of the user's action when s/he presses
                %"a","v" or "r". 
                parameters = [];
                %CG: parameters contains the data from each event - the amplitude, rise
                %time constant and decay time constant. These parameters are arranged in
                %columns and align with the single string values "a","v" or "r". Data is
                %ordered in a single column as follows:
                %parameters(1,x) = ROI index value (find roi in CA_ROIs etc)
                %parameters(2,x) = amplitude value
                %parameters(3,x) = rise time constant [AA(3)]
                %parameters(4,x) = decay time constant [AA(2)]
                %parameters(5,x) = complete trace event detection baseline indices
                %parameters(6,x) = baseline detection y-values
                %parameters(7,x) = complete trace event detection peak indices
                %parameters(8,x) = peak detection y-values

            end
            if ~isempty(undostring)
                a_indices = strfind(undostring, 'a'); acceptedtraces = numel(a_indices);
                v_indices = strfind(undostring, 'v'); partacceptedtraces = numel(v_indices);
                r_indices = strfind(undostring, 'r'); rejectedtraces = numel(r_indices);
                reviewedtraces = size(undostring,2);
            end
            accept_text.String = strcat('"', num2str(acceptedtraces), '" traces accepted');
            acceptamponly_text.String = strcat('"', num2str(partacceptedtraces), '" traces accepted (except curve fit)');
            reject_text.String = strcat('"', num2str(rejectedtraces), '" traces rejected');

            for cROI = 1 : Sz_CA_Measure(1);

                if ~isempty(CA_Measure{cROI,11})

                    CA_ew = CA_Measure{cROI,11}; 
%CG: numpotseg means 'number of potential segements'. CA_ew is a cell array
%containing all potential segements for the current ROI. A segement is a
%small section of a trace containing a single event. 
                    numpotseg = size(CA_ew,2); seg_indices = [];
                    for cRow = 1 : numpotseg
                        if ~isempty(CA_ew{cRow})
                            totaltraces = totaltraces + 1; seg_indices = [seg_indices; cRow];
                            if ~isempty(undostring)
                                if numel(undostring)+1 == totaltraces
                                    roi_idx = numel(roi_indices)+1; %trace_idx = numel(seg_indices); 
                                    trace_idx = numel(seg_indices); reviewedtraces = totaltraces;
                                end
                            end
                        end
                    end
                    if ~isempty(seg_indices)
                        roi_indices = [roi_indices; cROI]; CA_Segs{cROI} = seg_indices;
                    end
                    if cROI == roi_indices(roi_idx); 
                        CA_first = CA_ew; 
                    end
                    if totaltraces == numel(undostring) && isempty(roi_idx) 
                        roi_idx = numel(roi_indices); CA_first = CA_Measure{end,11}; 
                        counter = 1;
                        while isempty(CA_first)
                            CA_first = CA_Measure{end-counter,11};
                            counter = counter + 1;
                        end
                    end
                end
            end
            if isempty(undostring)
                roi_idx = 1; CA_first = CA_Measure{roi_indices(roi_idx),11};
%CG: roi_indices(roi_idx) correctly indexes the required segement. In this
%condition, there is no previous work to reload, so we must start from the
%first non-empty cell in the list.
            end
            %CG: some cells in CA_Measure and CA_ewmore arrays can be empty, so the
            %arrays roi_indices and seg_indices are used to keep a record of where the
            %non-empties are. 
            %CG: display number of segments

%             fh.Visible = 'on';
            newtext_top = text_top-(textvertspace*3);

            count_text_loc = round([(figpos(1)*textleft),...
                figpos(2)*newtext_top,winwidth*textwidth,winheight*textheight]);

            count_text = uicontrol('Style','text',...
                'Position',count_text_loc, 'FontSize', fontsz,...
                'HorizontalAlignment', 'left');

            count_text.String = strcat('"', num2str(reviewedtraces),'/',num2str(totaltraces),'" segments reviewed');

            seg_indices = CA_Segs{roi_indices(roi_idx)}; numrois = numel(roi_indices); numsegs = numel(seg_indices);
            CA_ew = CA_first; AmpVals = CA_Measure{roi_indices(roi_idx),9};
            ROI_IdxVal = CA_Measure{roi_indices(roi_idx),2};
            current_trace = CA_first{seg_indices(trace_idx)};

            Y = current_trace; X=[1:numel(Y)]'; 
            [ii, jj]=max(Y); 
            kk=min(Y); ll=find(Y>mean(Y));ll=numel(ll); 
            A  = [ii-kk, ll, ll/10, kk, jj, 0];
            F=[1 1 1 1 1 0]; E=F*0; [chi, N]=FIT_K_2EXP_SHIFT(X, Y, A, F ,E);
            figure(fh); plot(X,Y); 
            hold on; AA=A; X=[1:0.1:X(end)]; 
            Yf = myEXP_Shift([AA(1) AA(2) 0 AA(5)],X)-myEXP_Shift([AA(1) AA(3) 0 AA(5)],X) + AA(4);
            plot(X,Yf); 
            %AA(2) = T2; decay time constant
            %AA(3) = T1; rise time constant

            %CG: plot markers for baseline and peak for the first trace
            ybases = CA_Measure{roi_indices(roi_idx),7}; ypeaks = CA_Measure{roi_indices(roi_idx),8};
            xbases = CA_Measure{roi_indices(roi_idx),3}; xpeaks = CA_Measure{roi_indices(roi_idx),4};

            xb = xbases(seg_indices(trace_idx)); yb = ybases(seg_indices(trace_idx));
            xp = xpeaks(seg_indices(trace_idx)); yp = ypeaks(seg_indices(trace_idx));

            plot(xb,yb,'.','Color', [0 0.75 0], 'MarkerSize', 18)
            plot(xp,yp,'r.','MarkerSize', 18); hold off
            CA_ewidx = CA_Measure{roi_indices(roi_idx),12}; x_realframe = CA_ewidx{seg_indices(trace_idx)};

            ce = 1; CA_x_realframe = {}; breakWhileLoop = 0;
            while breakWhileLoop == 0 
                CA_x_realframe{1,ce} = x_realframe(str2double(ha.XTickLabels{ce})+1);
                ce = ce + 1;
                if numel(ha.XTickLabels) < ce; breakWhileLoop = 1; end
                if breakWhileLoop == 0; 
                    if (numel(x_realframe) < str2double(ha.XTickLabels{ce})+1)
                        breakWhileLoop = 1;
                    end
                end
            end
            ha.XTickLabels = CA_x_realframe;
            ha.XTick(numel(CA_x_realframe)+1:end) = [];
            %CG: get the indices for each detection within the current roi (according
            %to the complete fluorescence signal - i.e. whole trace). Just for
            %verification purposes - see test_AXC function.
            xbases_cfs = CA_Measure{roi_indices(roi_idx),5}; xpeaks_cfs = CA_Measure{roi_indices(roi_idx),6};
            xb_cfs = xbases_cfs(seg_indices(trace_idx)); xp_cfs = xpeaks_cfs(seg_indices(trace_idx));

            movegui(fh,'center'); 
            ha.XLabel.String = 'number of frames'; ha.XLabel.FontSize = 14;
            ha.YLabel.String = 'raw fluorescence units'; ha.YLabel.FontSize = 14;



            %CG: if buttons have been turned off, the user can use keys to carry out
            %analyses.
            if totaltraces == numel(undostring);
                breakloop = 1;
            else
                breakloop = 0;
            end
            %breakloop = 0;
            while breakloop == 0 && cConditionIdx <= NumConds
%CG: cConditionIdx <= NumConds will prevent ReviewFit from attempting to
%load data from the next condition if there is one. 
                fh.Visible = 'on';
                try 
                    kk = waitforbuttonpress;
                catch 
                    text_status.String = 'Data review cancelled!';
                    breakloop = 1; kk = -1; cConditionIdx = NumConds + 1;
                end
                if kk == 1
                    KeyPressRes = get(fh,'CurrentCharacter');
                    if strcmp(KeyPressRes, 'a') || strcmp(KeyPressRes, 'r') || strcmp(KeyPressRes, 'v') || strcmp(KeyPressRes, 'M')

                        if (roi_idx == numrois && trace_idx == numsegs) || strcmp(KeyPressRes, 'M')
            %CG: all trace segments from each roi have been reviewed.
                            breakloop = 1;
                        end

                        if strcmp(KeyPressRes, 'a')
                            undostring = strcat(undostring,'a');
                            acceptedtraces = acceptedtraces+1;
                            accept_text.String = strcat('"', num2str(acceptedtraces), '" traces accepted');

                            parameters(1,end+1) = ROI_IdxVal;
                            AmpVals = CA_Measure{roi_indices(roi_idx),9};
                            parameters(2,end) = AmpVals(seg_indices(trace_idx));
                            parameters(3,end) = AA(3);
                            parameters(4,end) = AA(2)/frameRate;
                            parameters(5,end) = xb_cfs;
                            parameters(6,end) = yb;
                            parameters(7,end) = xp_cfs;
                            parameters(8,end) = yp;

                        elseif strcmp(KeyPressRes, 'r')
                            undostring = strcat(undostring,'r');
            %CG: the user hits the 'r' key. This is equivalent to clicking the 'Reject' button.
                            rejectedtraces = rejectedtraces+1;
                            reject_text.String = strcat('"', num2str(rejectedtraces), '" traces rejected');

                            parameters(1,end+1) = ROI_IdxVal;
                            parameters(2,end) = NaN;
                            parameters(3,end) = NaN;
                            parameters(4,end) = NaN;
                            parameters(5,end) = NaN;
                            parameters(6,end) = NaN;
                            parameters(7,end) = NaN;
                            parameters(8,end) = NaN;
                        elseif strcmp(KeyPressRes, 'v')
                            undostring = strcat(undostring,'v');
            %CG: accept amplitude measurement only. This means the user thinks the baseline (green marker)
            %and the peak (red marker) have been appropriately positioned. Equivalent
            %to clicking 'Accept Ampl. Only'.
                            partacceptedtraces = partacceptedtraces+1;
                            acceptamponly_text.String = strcat('"', num2str(partacceptedtraces), '" traces accepted (except curve fit)');

                            parameters(1,end+1) = ROI_IdxVal;
                            AmpVals = CA_Measure{roi_indices(roi_idx),9};
                            parameters(2,end) = AmpVals(seg_indices(trace_idx));
                            parameters(3,end) = NaN;
                            parameters(4,end) = NaN;
                            parameters(5,end) = xb_cfs;
                            parameters(6,end) = yb;
                            parameters(7,end) = xp_cfs;
                            parameters(8,end) = yp;
                        end
                        save(filestr, 'undostring', 'parameters','-append')
                        if breakloop == 0;

                            if roi_idx <= numrois && trace_idx < numsegs
            %CG: still traces to check in current roi. This means we stay on the current 
            %ROI but jump to the next trace segment.
                                trace_idx = trace_idx + 1;
                            elseif roi_idx < numrois && trace_idx == numsegs 
            %CG: all traces in the current roi have been reviewed. Must switch to next
            %roi and start with first trace in the list. 
                                roi_idx = roi_idx+1; trace_idx = 1;
                                ybases = CA_Measure{roi_indices(roi_idx),7}; ypeaks = CA_Measure{roi_indices(roi_idx),8};
                                xbases = CA_Measure{roi_indices(roi_idx),3}; xpeaks = CA_Measure{roi_indices(roi_idx),4};
                                xbases_cfs = CA_Measure{roi_indices(roi_idx),5}; xpeaks_cfs = CA_Measure{roi_indices(roi_idx),6};
                                CA_ew = CA_Measure{roi_indices(roi_idx),11};
                                seg_indices = CA_Segs{roi_indices(roi_idx)}; numsegs = numel(seg_indices);
                                ROI_IdxVal = CA_Measure{roi_indices(roi_idx),2};
                                CA_ewidx = CA_Measure{roi_indices(roi_idx),12};
                            end 
                            xb = xbases(seg_indices(trace_idx)); yb = ybases(seg_indices(trace_idx));
                            xp = xpeaks(seg_indices(trace_idx)); yp = ypeaks(seg_indices(trace_idx));
                            xb_cfs = xbases_cfs(seg_indices(trace_idx)); xp_cfs = xpeaks_cfs(seg_indices(trace_idx));
                            current_trace = CA_ew{seg_indices(trace_idx)};
                            x_realframe = CA_ewidx{seg_indices(trace_idx)};

                            Y = current_trace; X=[1:numel(Y)]'; 
                            [ii, jj]=max(Y); 
                            kk=min(Y); ll=find(Y>mean(Y));ll=numel(ll); 
                            A  = [ii-kk, ll, ll/10, kk, jj, 0];
                            F=[1 1 1 1 1 0]; 
                            E=F*0; 
                            [chi, N]=FIT_K_2EXP_SHIFT(X, Y, A, F ,E);
                            figure(fh); plot(X,Y); 
                            hold on;
                             
                            ha.XLabel.String = 'number of frames'; ha.XLabel.FontSize = 14;
                            ha.YLabel.String = 'raw fluorescence units'; ha.YLabel.FontSize = 14;
                            AA=A; X=[1:0.1:X(end)]; 
                            Yf = myEXP_Shift([AA(1) AA(2) 0 AA(5)],X)-myEXP_Shift([AA(1) AA(3) 0 AA(5)],X) + AA(4);
                            plot(X,Yf);

                            ce = 1; CA_x_realframe = {}; breakWhileLoop = 0;
                            while breakWhileLoop == 0 
                                CA_x_realframe{1,ce} = x_realframe(str2double(ha.XTickLabels{ce})+1);
                                ce = ce + 1;
                                if numel(ha.XTickLabels) < ce; breakWhileLoop = 1; end
                                if breakWhileLoop == 0; 
                                    if (numel(x_realframe) < str2double(ha.XTickLabels{ce})+1)
                                        breakWhileLoop = 1;
                                    end
                                end
                            end
                            ha.XTickLabels = CA_x_realframe;
                            ha.XTick(numel(CA_x_realframe)+1:end) = [];
%                             for ce = 1 : numel(ha.XTickLabels)
%                                 if str2double(xTickArray{ce}) == 0
%                                     CA_x_realframe{1,ce} = x_realframe(1);
%                                 elseif str2double(xTickArray{ce}) > 0 && str2double(xTickArray{ce}) <= numel(x_realframe) 
%                                     CA_x_realframe{1,ce} = x_realframe(str2double(xTickArray{ce}));
%                                 elseif str2double(xTickArray{ce}) > numel(x_realframe)
%                                     CA_x_realframe{1,ce} = [];
%                                 end
%                             end
%                             ha.XTickLabels = CA_x_realframe;

                            plot(xb,yb,'.','Color', [0 0.75 0], 'MarkerSize', 18)
                            plot(xp,yp,'r.','MarkerSize', 18); hold off

                            reviewedtraces = reviewedtraces + 1;
                            count_text.String = strcat('"', num2str(reviewedtraces),'/',num2str(totaltraces),'" segments reviewed');
                        end
                    elseif strcmp(KeyPressRes, 'U')
                        if ~isempty(undostring)
                            if strcmp(undostring(end),'a') || strcmp(undostring(end),'r') || strcmp(undostring(end),'v')
            %CG: last action was an 'a' key press. 
                                if roi_idx >= 1 && trace_idx > 1
                                    trace_idx = trace_idx-1; 

                                elseif roi_idx > 1 && trace_idx == 1
                                    roi_idx = roi_idx-1; CA_ew = CA_Measure{roi_indices(roi_idx),11};
                                    seg_indices = CA_Segs{roi_indices(roi_idx)}; numsegs = numel(seg_indices);
                                    trace_idx = numsegs; ROI_IdxVal = CA_Measure{roi_indices(roi_idx),2};
                                    CA_ewidx = CA_Measure{roi_indices(roi_idx),12};

                                    ybases = CA_Measure{roi_indices(roi_idx),7}; ypeaks = CA_Measure{roi_indices(roi_idx),8};
                                    xbases = CA_Measure{roi_indices(roi_idx),3}; xpeaks = CA_Measure{roi_indices(roi_idx),4};
                                    xbases_cfs = CA_Measure{roi_indices(roi_idx),5}; xpeaks_cfs = CA_Measure{roi_indices(roi_idx),6};
                                end

                                current_trace = CA_ew{seg_indices(trace_idx)}; x_realframe = CA_ewidx{seg_indices(trace_idx)};
                                xb = xbases(seg_indices(trace_idx)); yb = ybases(seg_indices(trace_idx));
                                xp = xpeaks(seg_indices(trace_idx)); yp = ypeaks(seg_indices(trace_idx));
                                xb_cfs = xbases_cfs(seg_indices(trace_idx)); xp_cfs = xpeaks_cfs(seg_indices(trace_idx));

                                Y = current_trace; X=[1:numel(Y)]'; 
                                [ii, jj]=max(Y); 
                                kk=min(Y); ll=find(Y>mean(Y));ll=numel(ll); 
                                A  = [ii-kk, ll, ll/10, kk, jj, 0];
                                F=[1 1 1 1 1 0]; E=F*0; [chi, N]=FIT_K_2EXP_SHIFT(X, Y, A, F ,E);
                                figure(fh); plot(X,Y); hold on;
                                
                                ha.XLabel.String = 'number of frames'; ha.XLabel.FontSize = 14;
                                ha.YLabel.String = 'raw fluorescence units'; ha.YLabel.FontSize = 14;
                                AA=A; X=[1:0.1:X(end)]; 
                                Yf = myEXP_Shift([AA(1) AA(2) 0 AA(5)],X)-myEXP_Shift([AA(1) AA(3) 0 AA(5)],X) + AA(4);
                                plot(X,Yf); 

                                ce = 1; CA_x_realframe = {}; breakWhileLoop = 0;
                                while breakWhileLoop == 0 
                                    CA_x_realframe{1,ce} = x_realframe(str2double(ha.XTickLabels{ce})+1);
                                    ce = ce + 1;
                                    if numel(ha.XTickLabels) < ce; breakWhileLoop = 1; end
                                    if breakWhileLoop == 0; 
                                        if (numel(x_realframe) < str2double(ha.XTickLabels{ce})+1)
                                            breakWhileLoop = 1;
                                        end
                                    end
                                end
                                ha.XTickLabels = CA_x_realframe;
                                ha.XTick(numel(CA_x_realframe)+1:end) = [];
%                                 for ce = 1 : numel(ha.XTickLabels)
%                                     if str2double(xTickArray{ce}) == 0
%                                         CA_x_realframe{1,ce} = x_realframe(1);
%                                     elseif str2double(xTickArray{ce}) > 0 && str2double(xTickArray{ce}) <= numel(x_realframe) 
%                                         CA_x_realframe{1,ce} = x_realframe(str2double(xTickArray{ce}));
%                                     elseif str2double(xTickArray{ce}) > numel(x_realframe)
%                                         CA_x_realframe{1,ce} = [];
%                                     end
%                                 end
%                                 for ce = 1 : numel(x_realframe)
%                                     CA_x_realframe{1,ce} = x_realframe(ce);
%                                 end
                                ha.XTickLabels = CA_x_realframe;

                                plot(xb,yb,'.','Color', [0 0.75 0], 'MarkerSize', 18)
                                plot(xp,yp,'r.','MarkerSize', 18); hold off

                                reviewedtraces = reviewedtraces - 1;
                                count_text.String = strcat('"', num2str(reviewedtraces),'/',num2str(totaltraces),'" segments reviewed');

                                if strcmp(undostring(end),'a')
                                    acceptedtraces = acceptedtraces-1;
                                    accept_text.String = strcat('"', num2str(acceptedtraces), '" traces accepted');
                                elseif strcmp(undostring(end),'r')
                                    rejectedtraces = rejectedtraces-1;
                                    reject_text.String = strcat('"', num2str(rejectedtraces), '" traces rejected');
                                elseif strcmp(undostring(end),'v')
                                    partacceptedtraces = partacceptedtraces-1;
                                    acceptamponly_text.String = strcat('"', num2str(partacceptedtraces), '" traces accepted (except curve fit)');
                                end
                                parameters(:,end) = [];
                                undostring(end) = '';
                            end
                        elseif isempty(undostring) 
                           disp('no more commands to undo')
                        end
                    elseif strcmp(KeyPressRes, 'R')
            %CG: if user wishes to restart.
                        choice = questdlg('Restarting cannot be undone. Are you sure you want to restart? ', ...
                            '','Restart', 'Cancel restart','Cancel restart');
                        % Handle response
                        switch choice
                            case 'Restart'
                                trace_idx = 1; roi_idx = 1; CA_ew = CA_Measure{roi_indices(roi_idx),11};
                                seg_indices = CA_Segs{roi_indices(roi_idx)}; numsegs = numel(seg_indices);
                                current_trace = CA_ew{seg_indices(trace_idx)};
                                undostring = ''; parameters = [];
                                CA_ewidx = CA_Measure{roi_indices(roi_idx),12};
                                x_realframe = CA_ewidx{seg_indices(trace_idx)};

                                ybases = CA_Measure{roi_indices(roi_idx),7}; ypeaks = CA_Measure{roi_indices(roi_idx),8};
                                xbases = CA_Measure{roi_indices(roi_idx),3}; xpeaks = CA_Measure{roi_indices(roi_idx),4};
                                xb = xbases(seg_indices(trace_idx)); yb = ybases(seg_indices(trace_idx));
                                xp = xpeaks(seg_indices(trace_idx)); yp = ypeaks(seg_indices(trace_idx));

                                AmpVals = CA_Measure{roi_indices(roi_idx),9};

                                Y = current_trace; X=[1:numel(Y)]'; 
                                [ii, jj]=max(Y); 
                                kk=min(Y); ll=find(Y>mean(Y));ll=numel(ll); 
                                A  = [ii-kk, ll, ll/10, kk, jj, 0];
            %                         A(1) = AmpVals(seg_indices(trace_idx));
                                F=[1 1 1 1 1 0]; E=F*0; [chi, N]=FIT_K_2EXP_SHIFT(X, Y, A, F ,E);
                                figure(fh); plot(X,Y); hold on; 
                                ha.XLabel.String = 'number of frames'; ha.XLabel.FontSize = 14;
                                ha.YLabel.String = 'raw fluorescence units'; ha.YLabel.FontSize = 14;
                                AA=A; X=[1:0.1:X(end)]; 
                                Yf = myEXP_Shift([AA(1) AA(2) 0 AA(5)],X)-myEXP_Shift([AA(1) AA(3) 0 AA(5)],X) + AA(4);
                                plot(X,Yf); 

                                ce = 1; CA_x_realframe = {}; breakWhileLoop = 0;
                                while breakWhileLoop == 0 
                                    CA_x_realframe{1,ce} = x_realframe(str2double(ha.XTickLabels{ce})+1);
                                    ce = ce + 1;
                                    if numel(ha.XTickLabels) < ce; breakWhileLoop = 1; end
                                    if breakWhileLoop == 0; 
                                        if (numel(x_realframe) < str2double(ha.XTickLabels{ce})+1)
                                            breakWhileLoop = 1;
                                        end
                                    end
                                end
                                ha.XTickLabels = CA_x_realframe;
                                ha.XTick(numel(CA_x_realframe)+1:end) = [];
%                                 for ce = 1 : numel(ha.XTickLabels)
%                                     if str2double(xTickArray{ce}) == 0
%                                         CA_x_realframe{1,ce} = x_realframe(1);
%                                     elseif str2double(xTickArray{ce}) > 0 && str2double(xTickArray{ce}) <= numel(x_realframe) 
%                                         CA_x_realframe{1,ce} = x_realframe(str2double(xTickArray{ce}));
%                                     elseif str2double(xTickArray{ce}) > numel(x_realframe)
%                                         CA_x_realframe{1,ce} = [];
%                                     end
%                                 end
%                                 for ce = 1 : numel(x_realframe)
%                                     CA_x_realframe{1,ce} = x_realframe(ce);
%                                 end
                                ha.XTickLabels = CA_x_realframe;

                                plot(xb,yb,'.','Color', [0 0.75 0], 'MarkerSize', 18)
                                plot(xp,yp,'r.','MarkerSize', 18); hold off

                                acceptedtraces = 0;
                                accept_text.String = strcat('"0" traces accepted');
                                partacceptedtraces = 0;
                                acceptamponly_text.String = strcat('"0" traces accepted (except curve fit)');
                                rejectedtraces = 0;
                                reject_text.String = strcat('"0" traces rejected');

                                reviewedtraces = 1;
                                count_text.String = strcat('"', num2str(reviewedtraces),'/',num2str(totaltraces),'" segments reviewed');
                            case 'Cancel Restart'                     
                        end

                    end
                end
            end
            cConditionIdx = cConditionIdx + 1;
        end
    end
    if kk ~= -1; close(fh); end
end



