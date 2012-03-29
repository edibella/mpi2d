function cinemriUpsampled=mpi_upsampleImages(cinemri1,upsampleFactor)

if(~exist('upsampleFactor'))
  upsampleFactor=2;
end

% added 6/29/04 evrd
[nrows ncols nFrames]= size(cinemri1);
%nrows=2*(2*nrows-1)-1;
%ncols=2*(2*ncols-1)-1; % if odd same? yes.
cinemriUpsampled=[];
%newImg=interp2(tmpImg,2,'cubic');  so size is 2*(2*nrows-1)-1 x
%keyboard
warning('off','MATLAB:intMathOverflow');
  for i=1:nFrames
%       cc=interp2(cinemri1(:,:,i),2,'cubic');
       cinemriUpsampled(:,:,i) = interp2(cinemri1(:,:,i),upsampleFactor-1,'cubic');
       %cinemriUpsampled(:,:,i) = cinemriUpsampled(:,:,i)-min(min(cinemriUpsampled(:,:,i)));
       %cinemriUpsampled(:,:,i) = cinemriUpsampled(:,:,i)/(max(max(cinemriUpsampled(:,:,i))));
       %figure(45);
       %imshow(cinemriUpsampled(:,:,i));
       %cinemriUpsampled(:,:,i) = medfilt2(cinemriUpsampled(:,:,i),[3 3]);
       %figure(46);
       %imshow(cinemriUpsampled(:,:,i));
       %cinemriUpsampled(:,:,i) = histeq(cinemriUpsampled(:,:,i));
       %figure(47);
       %imshow(cinemriUpsampled(:,:,i));
       %pause;
       %disp(['Done upsampling of frame ' int2str(i)]);
%       cinemriUpsampled = [cinemriUpsampled cc];
  end
% same as xi,yi=meshgrid(1:.25:60,...
%cinemriUpsampled=reshape(cinemriUpsampled,[nrows ncols nFrames]);

warning('on','MATLAB:intMathOverflow');
return
