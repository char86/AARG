function [CellNameList, medianVals] = LatShiftCheck(cDir,suffixStrs,...
    CellNameList,medianVals,figFate)

% Outline
% LatShiftCheck checks if shift measurement has already been done for the
% current experiment. If so, it prompts the user to specify whether or not
% the previous measurement should be overwritten. If no previous
% measurement exists or it should be overwritten, then GetLatShift is
% called which calls mySmooth2D_All3_fft to measure the level of shift (in
% pixels) along the x and y axis between the first frame of the first
% condition and the last frame of the last condition. 

% Author: Charlie J Gilbride
% version: 1.0.0

if nargin == 0
    datadirs = uipickfiles('prompt','Select data folder(s) to check lateral shift');
%  datadirs = {}; %datadirs{1} = '/Users/cjgilbride/Documents/DecJan2016_17/MATLAB/AARG_pub_datatest/CG0706173';
% datadirs{1} = '/Volumes/My Passport for Mac/DecJan2016_17/MATLAB/AARG_Data/ReAnalysis/TTAP2/CG0408172';
% datadirs{1} = '/Volumes/My Passport for Mac/DecJan2016_17/MATLAB/AARG_Data/ReAnalysis/TTAP2/CG0408173';
    
    suffixStrs = {'a_baseline','b_TTAP2'}; partnamesuff = '-file00';
%CG: In the case of large tiff files that have been split by the
%acquisition program, the user can define "partnamesuff" such that all
%segments of the split file are converted to .mat files.
elseif nargin > 0
    partnamesuff = '';
    datadirs = {}; datadirs{1} = cDir; cd(cDir); dirInfo = dir; 
    ItemList = {}; NumItems = size(dirInfo,1); itemcount = 0;
    for cItemIdx = 1 : NumItems
        cItem = dirInfo(cItemIdx).name; 
%CG: first get all tif files together.
        if ~isempty(strfind(cItem,'.tif'))
            itemcount = itemcount + 1; 
            ItemList{itemcount,1} = cItem;
        end
    end
%CG: isolate tif files with some unidentified characters between keyword
%and file format string. 
    NumTifs = size(ItemList,1); specialstrs = {}; specialcount = 0;
    for cTifIdx = 1 : NumTifs
        strcount = 0; strfound = 0;
        while  strfound == 0 && strcount < size(suffixStrs,2)
            strcount = strcount + 1;
            if isempty(strfind(ItemList{cTifIdx},strcat(suffixStrs{strcount},'.tif'))) &&...
                    ~isempty(strfind(ItemList{cTifIdx},suffixStrs{strcount}))
                    
                strfound = 1; 
            end
        end
        if strfound == 1; 
            specialcount = specialcount + 1; 
            specialstrs{specialcount,1} = ItemList{cTifIdx}; 
        end
    end
    NumSpecStrs = size(specialstrs,1); 
%CG: isolate the unidentified string and search from back of string to find
%first 0. Final strings should be identical once the first 0 is found. 
    IDcount = 0; unID = {}; 
    for cSpecIdx = 1 : NumSpecStrs
        strcount = 0; strfound = 0; cSpecStr = specialstrs{cSpecIdx};
        while strcount < size(suffixStrs,2) && strfound == 0
            strcount = strcount + 1;
            stridx = strfind(cSpecStr,strcat(suffixStrs{strcount}));
            if ~isempty(stridx)
                tifidx = strfind(cSpecStr,'.tif'); IDcount = IDcount + 1;
                strID = cSpecStr(stridx+numel(suffixStrs{strcount}):tifidx);
                zeroidx = strfind(strID,'0'); 
                if ~isempty(zeroidx)
                    strID = strID(1:zeroidx(end));
                    if IDcount == 1;
                        unID{IDcount} = strID;
                    else
                        if ~strcmp(unID{1},strID)
                            unID{IDcount} = strID;
                        end
                    end
                end
            end
        end
    end
%CG: unID will contain the string identifying split files.
    if ~isempty(unID)
        if size(unID,2) > 1
            errordlg('partnamesuff not found')
        elseif size(unID,2) == 1
            partnamesuff = unID{1};
        end
    end

end
screensize = get( groot, 'Screensize' ); framesel = 1; %CG: frame to be displayed. Default: 1.

if ~isempty(partnamesuff)
    NumExpr = strcat(partnamesuff, '\d*');
else 
    NumExpr = '';
end

NumConditions = size(suffixStrs,2);

if iscell(datadirs) && ~isempty(datadirs)
    NumDirs = size(datadirs,2); 
    for cDirIdx1 = 1 : NumDirs
        cd(datadirs{cDirIdx1}); [~, CellName, ~] = fileparts(datadirs{cDirIdx1});        
        ExpNameList = {}; numdirs = size(datadirs,2); CellNameList{end+1} = CellName;
%CG: create FileNameLists.mat file if not already existing.
        try
%CG: In theory it is good to save FileNameLists.mat because it is a
%light-weight file but it can take a bit of time for matlab to put it
%together. However, if FileNameLists is wrongly configured and the user
%runs AARGI again, this will throw an error.
            DoNotTry
%             lists = load(strcat(CellName, '_FileNameLists.mat'));
%             ExpNameList = lists.ExpNameList;
%             partnamenumList = lists.partnamenumList;
%             NumStrFilesList = lists.NumStrFilesList;
%             suffixStrs_copy = suffixStrs;
%             expCount = 1;
%             while ~isempty(suffixStrs_copy) && expCount <= size(ExpNameList,1)
%                 if ~isempty(strfind(ExpNameList{expCount},suffixStrs_copy{expCount}))
%                     suffixStrs_copy{expCount} = {[]};
%                 end
%                 expCount = expCount + 1;
%             end
%             if ~isempty(suffixStrs_copy)
%                 breakTry;
%             end
            
        catch
            CurrentDir = datadirs{cDirIdx1}; cd(CurrentDir);
            numExpNameListFiles = 0; dirInfo=dir; numfiles=size(dirInfo,1);                                                   

            for cFile = 1 : numfiles
                currentFile = dirInfo(cFile).name;
                CondFound = 0; cCond = 1;
                while cCond <= NumConditions && CondFound == 0
                    if ~isempty(strfind(currentFile, suffixStrs{cCond}))
                        CondFound = 1;
                    end
                    cCond = cCond + 1;
                end

                if ~isempty(strfind(currentFile, '.tif')) && CondFound == 1
                    if ~isempty(NumExpr)
                        partnameExpr = regexp(currentFile, NumExpr, 'match');
                    else
                        partnameExpr = '';
                    end

                    wlCount = 1; acceptFileName = 0;
%CG: all file names should have the condition string followed by '.tif' OR
%condition string followed by partnameExpr and then '.tif'.
                    while wlCount <= size(suffixStrs,2) && acceptFileName == 0
                        if ~isempty(strfind(currentFile, strcat(suffixStrs{wlCount},'.tif')))
                            acceptFileName = 1; 
                        elseif ~isempty(partnameExpr)
                            if ~isempty(strfind(currentFile, strcat(suffixStrs{wlCount},partnameExpr,'.tif')))
                                acceptFileName = 1;
                            end
                        end
                        wlCount = wlCount + 1;
                    end
                    if isempty(partnameExpr) && acceptFileName == 1

                        ExpName = currentFile;
                        numExpNameListFiles = numExpNameListFiles + 1;
                        ExpNameList{numExpNameListFiles, 1} = ExpName;

                        partnamenumList{numExpNameListFiles, 1} = 1;

                    elseif ~isempty(partnameExpr) && acceptFileName == 1

                        ExpName = currentFile;
                        numExpNameListFiles = numExpNameListFiles + 1;
                        ExpNameList{numExpNameListFiles, 1} = ExpName;

                        partnameStr = partnameExpr{1};
                        Nulls = strfind(partnameStr, '0');
                        if ~isempty(Nulls)
                            partnameStr(1:Nulls(end)) = [];
                        end
                        partnamenumList{numExpNameListFiles, 1} = str2double(partnameStr);

                    end
                end
            end
            ExpNameList_sorted = cell(numExpNameListFiles, 1);
            partnamenumList_sorted = cell(numExpNameListFiles, 1);
            dExpr1 = strcat('(?<=', partnamesuff, ')\d');

            NumStrFilesList = zeros(numExpNameListFiles,numdirs);

            cDirIdx3 = 1;
            for cItemIdx = 1 : numExpNameListFiles
                for cCell = 1 : NumConditions
                    fn = ExpNameList{cItemIdx,cDirIdx3};
                    if ~isempty(strfind(fn, suffixStrs{cCell}))
                        NumStrFilesList(cItemIdx,cDirIdx3) = cCell;
                    end
                end
            end

            cDirIdx3 = 1;
            for cCell = 1 : NumConditions
                CondIdx = find(NumStrFilesList(:,1) == cCell);
                LastCondIdx = CondIdx(end); NumCondIdx = numel(CondIdx);
                for cItemIdx = 1 : numExpNameListFiles

                    if ~isempty(strfind(ExpNameList{cItemIdx, cDirIdx3}, suffixStrs{cCell}))
                        cItem2sort_ExpName = ExpNameList{cItemIdx, cDirIdx3};
                        partnameExpr1 = regexp(cItem2sort_ExpName, dExpr1, 'match');
                        if isempty(partnameExpr1); partnameExpr1 = '1'; end

                        if ~isempty(partnameExpr1)
                            cPartNameNum = partnamenumList{cItemIdx, cDirIdx3};
                            SortIdx = (LastCondIdx-NumCondIdx)+cPartNameNum;
                            ExpNameList_sorted{SortIdx, cDirIdx3} = ExpNameList{cItemIdx, cDirIdx3};
                            partnamenumList_sorted{SortIdx, cDirIdx3} = partnamenumList{cItemIdx, cDirIdx3};
                        end
                    end
                end
            end
            ExpNameList = ExpNameList_sorted;
            partnamenumList = partnamenumList_sorted;
            save(strcat(CellName, '_FileNameLists.mat'),'ExpNameList','partnamenumList',...
                'numExpNameListFiles','NumStrFilesList')
        end

        file = ExpNameList{1,1};

        tifnameIdx = strfind(file,'.tif'); FileName = file(1:tifnameIdx-1);
        file = double(imread(file,framesel));
%CG: just in case imread imports all frames instead of just the one frame
%specified by framesel:
        if numel(size(file)) == 3; file = file(:,:,framesel); end
        fh = figure('Name', strcat('"',FileName, '"; frame "', num2str(framesel),'"'),'Visible','Off'); imagesc(file); colormap jet
        
        left = screensize(3)*0.2; bot = screensize(2);
        wid = screensize(3)*0.66; hei =  screensize(4);
        ax1 = gca; 
        set(ax1, 'box', 'off'); set(fh, 'OuterPosition', [left, bot, wid, hei]);
        set(fh, 'NumberTitle', 'Off')
        axis image; 

        try
            sk = load(strcat(CellName, '_shiftkey.mat'));
            shiftkey = sk.shiftkey;
            answer = questdlg('Overwrite existing shiftkey variable?', ...
                'Overwrite', ...
                'Yes','No','No');
            % Handle response
            switch answer
                case 'Yes'
                    owctrl = 1;
                case 'No'
                    owctrl = 2;
            end
            if owctrl == 1
                fh.Visible = 'On';
                
                shiftkey = getrect(ax1); 
                rectangle('Position',shiftkey, 'EdgeColor', 'w')
                save(strcat(CellName, '_shiftkey.mat'), 'shiftkey');
            end
            
        catch
            fh.Visible = 'On';
            shiftkey = getrect(ax1); rectangle('Position',shiftkey, 'EdgeColor', 'w')
            save(strcat(CellName, '_shiftkey.mat'), 'shiftkey');
        end
        [medianVals] = GetLatShift(CellName,suffixStrs,ExpNameList,...
            NumStrFilesList,framesel,FileName,medianVals,figFate);
        
    end
    if figFate == 1; close(fh); end
end

end

function [medianVals] = GetLatShift(CellName,suffixStrs,ExpNameList,...
    NumStrFilesList,framesel,FileName,medianVals,figFate)

sk = load(strcat(CellName, '_shiftkey.mat')); shiftkey = sk.shiftkey;
NumConds = size(suffixStrs,2); ShiftCheckArray = zeros(NumConds,2);

%CG: Aim is to check for shift between the first and last frame either
%across all conditions or within each condition. 
for cCondIdx = 1 : NumConds
    idx = find(NumStrFilesList==cCondIdx);
    if cCondIdx == 1
        ShiftCheckArray(1,1) = 1;
        ShiftCheckArray(1,2) = idx(end);
    elseif cCondIdx > 1
        ShiftCheckArray(cCondIdx,1) = idx(1);
        ShiftCheckArray(end,2) = idx(end);
    end
end
%CG: assuming the user wants acceptable lateral shift between the first
%frame of the first condition and the last frame of the last condition. 
startfile = ExpNameList{ShiftCheckArray(1,1)}; idx = 0;
endfile = ExpNameList{ShiftCheckArray(end,2)}; 

ff = [];
for ci = 1 : 3
    ii = double(imread(startfile,'Index',ci)); 
    if ci == 1; ff = ii; else ff = cat(3,ff,ii); end
end
lf = double(imread(endfile(:,:,end)));
shiftkey = round(shiftkey);
Ref = ff(shiftkey(2):shiftkey(2)+shiftkey(4),shiftkey(1):shiftkey(1)+shiftkey(3),:);
Ref = mean(Ref-mySmooth2D_All3_fft(Ref,7,7),3);
fh1 = figure('Name', strcat('"',FileName, '"; frame "', num2str(framesel),'"')); imagesc(Ref); colormap jet

In = lf(shiftkey(2):shiftkey(2)+shiftkey(4),shiftkey(1):shiftkey(1)+shiftkey(3));
In = In-mySmooth2D_All3_fft(In,7,7); Out = myShiftDetect2D(In, Ref); 

% plot(Out(:,1),'k')
medianvalx = median(Out(:,1));
% hold on
% plot(Out(:,2))
medianvaly = median(Out(:,2));
% hold off

medianVals = [medianVals;medianvalx, medianvaly];
save(strcat(CellName, '_shiftkey.mat'),'medianVals','-append')
if figFate == 1; pause(3); close(fh1); end
end