% Extract histogram at point I from integral histogram
%
% ih = iHist(H, rbr, cbr)
% ih = iHist(H, rbr, cbr, rtl, ctl) 
% ih = iHist(H, rbr, cbr, rtl, ctl, b);
%
% Extracts the histogram of an image region from an integral histogram.
% 
% In the simplest form the histogram of the region with top left corner
% (0,0) and bottom right corner (rbr, cbr) is returned.
% 
% If b is specified, only the value of bin b is returned.
%
% If rtl, ctl is specified, it is used as the top left corner of the
% region.
%
% Ref: F. Porikli, Integral Histogram: a fast way to extract histograms in
% cartesian spaces".  Proc. Computer Vision an Pattern Recognition, vol 1,
% 829--836, 2008.
%

function ih = iHist(H, rbr, cbr, rtl, ctl, b) 

    if ~isa(H, 'cell')
        error('H must be a cell array of signed 32-bit histograms.');
    end

	% Bounds checking

    if rbr > size(H, 1)
        rbr = size(H,1);
    end
    
    if cbr > size(H,2)
        cbr = size(H,2);
    end
    
    if rtl > size(H,1)
        rtl = size(H,1);
    end
    
    if ctl > size(H, 2)
        ctl = size(H,2)
    end
    
    
    if nargin == 3

         ih = H{rbr,cbr};
       
    elseif nargin == 5 
        ih = H{rbr, cbr};
             
        if rtl > 1
            ih = ih - H{rtl-1, cbr};
        end
        
        if ctl > 1
            ih = ih - H{rbr, ctl-1};
        end

        if ctl > 1 && rtl > 1
            ih = ih + H{rtl-1, ctl-1};
        end
      
        
    elseif nargin == 6 
        
        ih = H{rbr,cbr}(b);
        
        if rtl > 1
            ih = ih - H{rtl-1,cbr}(b);
        end
        
        if ctl > 1
            ih = ih - H{rbr,ctl-1}(b);
        end
        
        if ctl > 1 && rtl > 1
            ih = ih + H{rtl-1, ctl-1}(b);
        end
     
        
    else

        error('Incorrect number of arguments to iHist.');
        
    end

    
end