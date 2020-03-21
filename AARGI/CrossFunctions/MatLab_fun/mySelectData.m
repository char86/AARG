function varargout = mySelectData(nPos, nMore, Str1)
% Out = mySelectData(nPos, nMore, Str);
% e.g.
% Out = mySelectData([], 1);
% or
% Out = mySelectData([], 1, Str); % Str(:).name;  Str(:).size; Str(:).class  will be displaied

if nargin<3; ntt = 1; wsVar = evalin('base','whos'); else wsVar = Str1; ntt=2; end
gl_ii=[];
while 1
    if numel(wsVar)<1; if ntt==1; warndlg('The workspace is empty !'); else warndlg('The choose-item is empty !');  end; break; end
    %sz=zeros(numel(wsVar),1); for ni=1:numel(wsVar); sz(ni) = numel(wsVar(ni).name); end; mm=max(sz(:));
    if ntt==1
        for ni=1:numel(wsVar);
            %str{ni} = ['<html>', wsVar(ni).name, repmat('&nbsp',[1  mm-sz(ni)]), '<FONT color="#0000ff"><i>',...
            %    mySize2Str_ex(wsVar(ni).size),' ',wsVar(ni).class ,'</FONT>'];
            str{ni} = ['<html><table><tr><td width=77>', wsVar(ni).name, '</td><td><FONT color="#176797"><i>',...
                mySize2Str_ex(wsVar(ni).size),' ',wsVar(ni).class ]; %mySize2Str_ex(wsVar(ni).size),' ',wsVar(ni).class ,'</i></FONT></td></table></htm;>'];
        end %for ni=1:numel(wsVar); str{ni} = wsVar(ni).name; end
    else
        for ni=1:numel(wsVar);
            str{ni} = ['<html><table><tr><td width=77>', wsVar(ni).name, '</td><td><FONT color="#176797"><i>',...
                mySize2Str_ex(wsVar(ni).size),' ',wsVar(ni).class ]; %mySize2Str_ex(wsVar(ni).size),' ',wsVar(ni).class ,'</i></FONT></td></table></htm;>'];
        end %for ni=1:numel(wsVar); str{ni} = wsVar(ni).name; end
    end
    %if nargin<1; nPos = [100,100,200,270]; elseif numel(nPos)<1; nPos = [100,100,200,270]; end
    if nargin<1; nPos = myDefaultPos(ntt, wsVar, 1); elseif numel(nPos)<1; nPos = myDefaultPos(ntt, wsVar, 2); end
    if nargin<2; nMore=0; elseif numel(nMore)<1; nMore=0; end
    gl_ss = get(0,'screensize'); nPos(2)=gl_ss(4)-nPos(4)-nPos(2); % Screen -> matlab position
    hh = figure('Position',nPos,'MenuBar','none','Name','Select a Valable',...
        'NumberTitle','off','Resize','off','CloseRequestFcn',@myCallBack_Close);
    nP=[10,50,nPos(3)-20,nPos(4)-60];
    hList = uicontrol('Style','listbox','Position', nP, 'String',str,'Value',1); nP=[10,10,nPos(3)/2-15,20];
    if nMore>0; hList.Max=hList.Min+2;end
    uicontrol('Style','pushbutton','Position',nP,'String','OK','FontWeight','bold','Callback',@myCallBack_OK);
    nP=[5+nPos(3)/2,10,nPos(3)/2-15,20];
    uicontrol('Style','pushbutton','Position',nP,'String','Cancel', 'FontWeight','bold','Callback',@myCallBack_Close);
    try
        uicontrol(hList); uiwait(hh);
    catch
        if ishghandle(hh); delete(hh); end
    end
    break;
end
%if nargout>0; if numel(gl_ii)>0; varargout{1} = evalin('base',myCombine_ex(str, gl_ii)); else  varargout{1} =[]; end; end
if ntt==1
    if nargout>0; if numel(gl_ii)>0; varargout{1} = evalin('base',myCombine_ex(wsVar, gl_ii)); else  varargout{1} =[]; end; end
else
    if nargout>0;
        if numel(gl_ii)>0;
            if numel(gl_ii)>1;
                ni = gl_ii(1); nj=ni;
                for nk=2:numel(gl_ii)
                    if numel(wsVar(ni).size)==numel(wsVar(gl_ii(nk)).size)
                        if sum(wsVar(ni).size==wsVar(gl_ii(nk)).size)==numel(wsVar(ni).size)
                            nj=[nj gl_ii(nk)];
                        end
                    end
                end
                gl_ii=nj;
            end
            varargout{1} = gl_ii;
        else  varargout{1} =[]; end; end
end
drawnow; % Update the view to remove the closed figure (g1031998)

    function myCallBack_OK(varargin)
        gl_ii = hList.Value;
        delete(gcbf);
    end
    function myCallBack_Close(varargin); delete(gcbf); end
    function Out = myCombine_ex_o(In, jj)
        Out='['; for ii=1:numel(jj); Out=[Out ' ' In{jj(ii)}]; end; Out=[Out ' ]'];
    end
    function Out = myCombine_ex(In, jj)
        Out='['; for ii=1:numel(jj); Out=[Out ' ' In(jj(ii)).name]; end; Out=[Out ' ]'];
    end
    function str = mySize2Str_ex(nn)
        if numel(nn)<1;str=''; else
            str=num2str(nn(1)); for ii=2:numel(nn); str=cat(2, str,['x', num2str(nn(ii))] ); end
        end
    end
    function Out=myDefaultPos(ntt, wsVar, pp)
        if nargin<3; pp=1; end
        if ntt==1; if pp==1; Out = [100,100,200,270]; else  Out = [100,100,200,270]; end
        else aa = numel([wsVar(1).name wsVar(1).class])*8 + (4*4+7)*7;
            Out = [100,100,200,270];
            if(Out(3)<aa); Out(3) = aa; end
        end
    end
end
