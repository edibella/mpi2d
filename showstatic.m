%
% Compute static shift for 
% the whole image
%
function showstatic(m,m2,delay, referenceFrame)

[nx,ny,nt]=size(m);

for it=1:nt-1

 r1=(m(:,:,it+1)-m(:,:,it));
 r2=(m2(:,:,it+1)-m2(:,:,it));
 r3=(m2(:,:,it)-m(:,:,referenceFrame));

 %subplot(2,2,1); imagesc(m(:,:,it)); title('Previos'); 
 %subplot(2,2,2); imagesc(m(:,:,it+1)); title('Next');
 subplot(2,2,1); imagesc(m(:,:,referenceFrame)); title('Reference frame,  no shifts '); 
 fprintf('Reference frame used is %d \n',referenceFrame); 
 subplot(2,2,2); imagesc(m2(:,:,it)); title('After shift');
 %subplot(2,2,3); imagesc(r1); title('diff consecutive Originals'); caxis(300*[-1 1]);
 subplot(2,2,3); imagesc(r2); title('diff consecutive Outputs'); caxis(300*[-1 1]);
 subplot(2,2,4); imagesc(r3); title('diff with Reference'); caxis(300*[-1 1]);

 disp([ it var( r1(:)) var( r2(:)) ]);
 %pause(delay);
 pause;

end

return;
