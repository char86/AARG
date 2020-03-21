function [Out mm]= myRemove_BK(In)
% Using:
% Out = myRemove_BK(In);
% or
% [Out, mm] = myRemove_BK(In);
%
% In is 3D, 2D-t, data is in last-DIM.
% Remove a global background in last-DIM
%
% On return:
% Out is BK removed data.
% mm is the signal removed.

ss = size(In);In2 = reshape(In,[],ss(end));
mm = mean(In2,1);
Out = reshape(myRemove_in(In2, mm),ss);

    function Out = myRemove_in(In, mm)
        % Out = myRemove(In, S)
        % In is 2d, last dim is time
        % S single line with length = last dim of In
        % remove the co-relation part of each pixel in In,
        mm = mm(:); 
%CG: mm = mean of all 28900 pixels per frame.        
        nn = numel(mm) - 1;
        mm = mm - mean(mm);
        mm = mm./sqrt(sum(mm.*mm)./nn);
%CG: now mm represents normalized average for each frame (dividing by
%std). mm(t) will be larger when more signal occurs in the current frame t.
        mx = mean(In, 2);
%CG: mx = mean of each pixel across time. mx(i) will be larger when the
%current pixel overlaps with an active part of the cell. 
        ddx = repmat(mx,1,nn + 1);
        X = In - ddx;
%CG: X contains the background of In 
        c = X * mm./nn;
%CG: c is the pixelwise cross-correlation coefficient between the input
%data with mean 0 (X) and the total average with mean 0 (mm). 
%CG: c(i) will be larger when background is higher (why divide by nn?). c
%is the correlated part of each pixel of In
        Out = In - c * mm';   
%CG: why multiply by mm again? 

    end
end