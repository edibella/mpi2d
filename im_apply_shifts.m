    function [imsr] =im_apply_shifts(ims, sh,flag);

% Image Registration: Apply Shifts to Register
%   
%   INPUT:
%       ims - A stack of images to be registred
%       sh  - Shifts in respect to the first image
%             
%   OUTPUT: 
%       imsr - Registred stack of images
%
%   by Alexey Samsonov
%   University of Wisconsin, Madison
%   January 2005

dims=size(ims);
if(flag==1)
nim=dims(3);
else
nim=1;
end

rx=linspace(-0.5, 0.5, dims(1)+1);
ry=linspace(-0.5, 0.5, dims(2)+1);
[xx yy]=ndgrid(rx(1:end-1),ry(1:end-1));

% correcting
imsr=ims;
for ii=1:nim
    ph=fftshift(-2*pi*(xx*sh(ii,1)+yy*sh(ii,2)));
    imsr(:,:,ii)=ifft2(fft2(squeeze(ims(:,:,ii))).*exp(i*ph));
end

return
