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
%CG: mm = mean of all 28900 pixels per frame. Similar to y~(i) in Chen et al (2006).        
        nn = numel(mm) - 1;
        mm = mm - mean(mm);
        mm = mm./sqrt(sum(mm.*mm)./nn);
%CG: now mm represents normalized average for each frame (dividing by
%std). mm is the same as the first component dervied from singular value
%decomposition. We use the mean here instead of svd (as in Chen et al
%(2006)) because the svd plus regression method is slower when there is
%only one ROI. In this analysis, we treat the entire field of view as a
%single ROI for which the background signal is removed. Chen et al (2006)
%in contrast, calculated the background signal for each ROI individually.
%Finding the background signal using the regression method would be slower
%for our data (i.e. single cells expressing GCaMP6s and imaged at high mag
%to measure spine calcium signals). Using the regression method is also not
%necessary because the background signal is not as inhomogenous as in Chen
%et al (2006). Chen et al (2006) did not use a GECI but rather a
%bath-applied membrane permeable dye to label lots of cells, so their
%background is higher and more inhomogenous. 
        mx = mean(In, 2);
%CG: mx = mean of each pixel across time. mx(i) will be larger when the
%current pixel overlaps with an active part of the cell. mx is an offset
%that is removed from the pixel signal. 
        ddx = repmat(mx,1,nn + 1);
        X = In - ddx;
%CG: X contains the background of In 
        c = X * mm./nn;
%CG: c is the pixelwise cross-correlation coefficient between the input
%data with mean 0 (X) and the total average with mean 0 (mm). 
%CG: c(i) will be larger when background is higher. The product of c*mm'
%will extract the noise component from the signal which can then be
%subtracted from the data matrix. 
        Out = In - c * mm';   

    end
end

%Reference: Chen et al (2006). PMID:16387783