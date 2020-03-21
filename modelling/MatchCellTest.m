function MatchCellTest()

TgtDir = uigetdir('','Select folder containing model data and AARG result');
% TgtDir = '/Users/cjgilbride/Documents/DecJan2016_17/MATLAB/AARG_pub/Latest/Modelling/5x5_ROIs/ModelDataFour_KR';
addpath(TgtDir);
cd(TgtDir)
[~, ModelDataDir, ~] = fileparts(TgtDir);
DefString = ModelDataDir;
load(strcat('SLM_CellArrays_', DefString, '.mat'), 'SLM_CA_ROIs', 'SLM_CA_EventNum',...
    'SLM_CA_EventSize', 'SLM_CA_FrameIEI', 'SLM_CA_EventType')

if strcmp(DefString(end-1:end), 'KR')
    load(strcat('PE_SLM_CAs_', DefString, '.mat'), 'SLM_CA_ROIs_PreEst', 'SLM_CA_EventSize_PreEst')
end

disp(strcat('MatchCellTest called for "', ModelDataDir, '"'));
cd(strcat(TgtDir, '/AARG_result'));
dirInfo = dir;                                                                
NumItems = size(dirInfo,1);
for cItemIdx = 1 : NumItems 
    cItem = dirInfo(cItemIdx).name;
    load(strcat('CAs_Core_', DefString, '.mat'), 'CA_ROIs', 'CA_NumEvents',...
       'CA_EventSize', 'CA_FrameIEI')
%     CtrlSwitch = 'Go!';
    if strfind(cItem, '_Core_')
       load(strcat('CAs_Core_', DefString, '.mat'), 'CA_ROIs', 'CA_NumEvents',...
           'CA_EventSize', 'CA_FrameIEI')
       disp('AARG_Core results...')
       CtrlSwitch = 'Go!';
    else
        CtrlSwitch = '';
    end

    if ~isempty(CtrlSwitch)
        Sz_SLM_CA_ROIs = size(SLM_CA_ROIs);
        Sz_CA_ROIs = size(CA_ROIs);

        if numel(Sz_SLM_CA_ROIs) ~= numel(Sz_CA_ROIs)
            errordlg('Modelling script and AARG do not agree on the dimensions of the CA_ROIs'); 
        end

        PDCounter_ROIs = 0;
        if numel(Sz_SLM_CA_ROIs) == 2   
            TotalCellNum_ROI = Sz_SLM_CA_ROIs(1)*Sz_SLM_CA_ROIs(2);
        elseif numel(Sz_SLM_CA_ROIs) == 3
            TotalCellNum_ROI = Sz_SLM_CA_ROIs(1)*Sz_SLM_CA_ROIs(2)*Sz_SLM_CA_ROIs(3);
        end
        for CurrentCell = 1 : TotalCellNum_ROI
            
            if ~isempty(SLM_CA_ROIs{CurrentCell}) && ~isempty(CA_ROIs{CurrentCell})
                R = SLM_CA_ROIs{CurrentCell} - CA_ROIs{CurrentCell};
                if R == 0
                    PDCounter_ROIs = PDCounter_ROIs + 1;
                elseif isempty(R)
                    PDCounter_ROIs = PDCounter_ROIs + 1;
                elseif R ~= 0
            %         CurrentCell
                    [Row_ROIs, Column_ROIs, Zspace_ROIs] = ind2sub(size(zeros(Sz_SLM_CA_ROIs)), CurrentCell);
                    disp(strcat('ROI error at: Frame ', num2str(Row_ROIs), ', component ',...
                        num2str(Column_ROIs), ', ZSpace ', num2str(Zspace_ROIs)));
                end
            elseif ~isempty(SLM_CA_ROIs{CurrentCell}) && isempty(CA_ROIs{CurrentCell}) || isempty(SLM_CA_ROIs{CurrentCell}) && ~isempty(CA_ROIs{CurrentCell})
                [Row_ROIs, Column_ROIs, Zspace_ROIs] = ind2sub(size(zeros(Sz_SLM_CA_ROIs)), CurrentCell);
                disp(strcat('ROI error at: Frame ', num2str(Row_ROIs), ', component ',...
                    num2str(Column_ROIs), ', ZSpace ', num2str(Zspace_ROIs)));
            elseif isempty(SLM_CA_ROIs{CurrentCell}) && isempty(CA_ROIs{CurrentCell}) 
                PDCounter_ROIs = PDCounter_ROIs + 1;
            end
        end

        disp(strcat('% agreement for ROI cell arrays :', num2str((PDCounter_ROIs/TotalCellNum_ROI)*100)))

        Sz_SLM_CA_EventNum = size(SLM_CA_EventNum);
        Sz_CA_NumEvents = size(CA_NumEvents);

        if numel(Sz_SLM_CA_EventNum) ~= numel(Sz_CA_NumEvents)
            errordlg('Modelling script and AARG do not agree on the dimensions of the CA_EventNum'); 
        end

        PDCounter_EN = 0;
        if numel(Sz_SLM_CA_EventNum) == 2   
            TotalCellNum_EN = Sz_SLM_CA_EventNum(1)*Sz_SLM_CA_EventNum(2);
        elseif numel(Sz_SLM_CA_EventNum) == 3
            TotalCellNum_EN = Sz_SLM_CA_EventNum(1)*Sz_SLM_CA_EventNum(2)*Sz_SLM_CA_EventNum(3);
        end
        for CurrentCell = 1 : TotalCellNum_EN

            if ~isempty(CA_NumEvents{CurrentCell}) && ~isempty(SLM_CA_EventNum{CurrentCell})
                R = SLM_CA_EventNum{CurrentCell} - CA_NumEvents{CurrentCell};

                if R == 0
                    PDCounter_EN = PDCounter_EN + 1;
                elseif isempty(R)
                    PDCounter_EN = PDCounter_EN + 1;
                elseif R ~= 0
                    [Row_Num, Column_Num, Zspace_Num] = ind2sub(size(zeros(Sz_SLM_CA_EventNum)), CurrentCell);
                    disp(strcat('Num error at: Frame ', num2str(Row_Num), ', component ',...
                        num2str(Column_Num), ', ZSpace ', num2str(Zspace_Num)));               
                end
            elseif ~isempty(CA_NumEvents{CurrentCell}) && isempty(SLM_CA_EventNum{CurrentCell}) || isempty(CA_NumEvents{CurrentCell}) && ~isempty(SLM_CA_EventNum{CurrentCell})
                [Row_Num, Column_Num, Zspace_Num] = ind2sub(size(zeros(Sz_SLM_CA_EventNum)), CurrentCell);
                disp(strcat('Num error at: Frame ', num2str(Row_Num), ', component ',...
                    num2str(Column_Num), ', ZSpace ', num2str(Zspace_Num)));
            elseif isempty(CA_NumEvents{CurrentCell}) && isempty(SLM_CA_EventNum{CurrentCell})
                PDCounter_EN = PDCounter_EN + 1;
            end
        end

        disp(strcat('% agreement for Event Number cell arrays :', num2str((PDCounter_EN/TotalCellNum_EN)*100)))

        Sz_SLM_CA_EventSize = size(SLM_CA_EventSize);
        Sz_CA_EventSize = size(CA_EventSize);

        if numel(Sz_SLM_CA_EventSize) ~= numel(Sz_CA_EventSize)
            errordlg('Modelling script and AARG do not agree on the dimensions of the CA_EventSize'); 
        end

        PDCounter_ES = 0;
        if numel(Sz_SLM_CA_EventSize) == 2   
            TotalCellNum_ES = Sz_SLM_CA_EventSize(1)*Sz_SLM_CA_EventSize(2);
        elseif numel(Sz_SLM_CA_EventSize) == 3
            TotalCellNum_ES = Sz_SLM_CA_EventSize(1)*Sz_SLM_CA_EventSize(2)*Sz_SLM_CA_EventSize(3);
        end
        for CurrentCell = 1 : TotalCellNum_ES

            if ~isempty(SLM_CA_EventSize{CurrentCell}) && ~isempty(CA_EventSize{CurrentCell})
                if numel(SLM_CA_EventSize{CurrentCell}) == numel(CA_EventSize{CurrentCell})
                    SLMans = SLM_CA_EventSize{CurrentCell}; AARGans = CA_EventSize{CurrentCell};
                    Sz_AARGans = size(AARGans);
                    if Sz_AARGans(1) > Sz_AARGans(2); AARGans = AARGans'; end
                    R = SLMans - AARGans;

                    if R == 0
                        PDCounter_ES = PDCounter_ES + 1;
                    elseif isempty(R)
                        PDCounter_ES = PDCounter_ES + 1;
                    end
                else
                    [Row_ES, Column_ES, Zspace_ES] = ind2sub(size(zeros(Sz_SLM_CA_EventSize)), CurrentCell);
                    disp(strcat('Size error at: Frame ', num2str(Row_ES), ', component ',...
                        num2str(Column_ES), ', ZSpace ', num2str(Zspace_ES)));
                end
            elseif isempty(SLM_CA_EventSize{CurrentCell}) && ~isempty(CA_EventSize{CurrentCell})
                PDCounter_ES = PDCounter_ES + 1;
            elseif isempty(SLM_CA_EventSize{CurrentCell}) && isempty(CA_EventSize{CurrentCell})
                PDCounter_ES = PDCounter_ES + 1;
            end

        end

        disp(strcat('% agreement for Event Size cell arrays :', num2str((PDCounter_ES/TotalCellNum_ES)*100)))

        Sz_SLM_CA_FrameIEI = size(SLM_CA_FrameIEI);
        Sz_CA_FrameIEI = size(CA_FrameIEI);

        if numel(Sz_SLM_CA_FrameIEI) ~= numel(Sz_CA_FrameIEI)
            errordlg('Modelling script and AARG do not agree on the dimensions of the CA_FrameIEI'); 
        end

        PDCounter_IEI = 0;
        if numel(Sz_SLM_CA_FrameIEI) == 2   
            TotalCellNum_IEI = Sz_SLM_CA_FrameIEI(1)*Sz_SLM_CA_FrameIEI(2);
        elseif numel(Sz_SLM_CA_FrameIEI) == 3
            TotalCellNum_IEI = Sz_SLM_CA_FrameIEI(1)*Sz_SLM_CA_FrameIEI(2)*Sz_SLM_CA_FrameIEI(3);
        end
        for CurrentCell = 1 : TotalCellNum_IEI
            
            if ~isempty(SLM_CA_FrameIEI{CurrentCell}) && ~isempty(CA_FrameIEI{CurrentCell})
                R = SLM_CA_FrameIEI{CurrentCell} - CA_FrameIEI{CurrentCell};

                if R == 0
                    PDCounter_IEI = PDCounter_IEI + 1;
                elseif isempty(R)
                    PDCounter_IEI = PDCounter_IEI + 1;
                elseif R ~= 0
                    [Row_IEI, Col_IEI, ZSpace_IEI] = ind2sub(size(zeros(Sz_SLM_CA_FrameIEI)), CurrentCell);
                    disp(strcat('IEI error at: Frame ', num2str(Row_IEI), ', component ',...
                        num2str(Col_IEI), ', ZSpace ', num2str(ZSpace_IEI)));
                end
            elseif ~isempty(SLM_CA_FrameIEI{CurrentCell}) && isempty(CA_FrameIEI{CurrentCell}) || isempty(SLM_CA_FrameIEI{CurrentCell}) && ~isempty(CA_FrameIEI{CurrentCell})
                [Row_IEI, Col_IEI, ZSpace_IEI] = ind2sub(size(zeros(Sz_SLM_CA_FrameIEI)), CurrentCell);
                disp(strcat('IEI error at: Frame ', num2str(Row_IEI), ', component ',...
                    num2str(Col_IEI), ', ZSpace ', num2str(ZSpace_IEI)));
            elseif isempty(SLM_CA_FrameIEI{CurrentCell}) && isempty(CA_FrameIEI{CurrentCell})
                PDCounter_IEI = PDCounter_IEI + 1;
            end

        end

        disp(strcat('% agreement for IEI cell arrays :', num2str((PDCounter_IEI/TotalCellNum_IEI)*100)))
    end
end    




