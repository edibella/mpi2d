% Integral histogram computation.
%
% Computes the integral histogram of 8-bit grayscale image I.
%
% Ref: F. Porikli, Integral Histogram: a fast way to extract histograms in
% cartesian spaces".  Proc. Computer Vision an Pattern Recognition, vol 1,
% 829--836, 2005.
%




function H = imIntegralHist(I)

    [height, width] = size(I);

    if ~isa(I, 'uint8')
        disp('imIntegralHist: Input image must be 8-bit grayscale');
    end

    H = cell(height, width);
      
    for r = 1:height
        
        for c = 1:width

              
              H{r,c} = zeros(1,256, 'int32');
              
             
              if r > 1 
                  H{r,c} = H{r,c} + H{r-1,c};                  
              end
              if c > 1
                  H{r,c} = H{r,c} + H{r, c-1};                  
              end

              if c > 1 && r > 1
                  H{r,c} = H{r,c} - H{r-1, c-1};
              end
              
              H{r,c}(floor(I(r,c)+1)) = H{r,c}(floor(I(r,c)+1)) + 1;
              
        end
    end
       
    
    