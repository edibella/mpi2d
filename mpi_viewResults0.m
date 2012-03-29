function mpi_viewResults(curves, studyNum, numSlices, outpath, numAzimuthalRegions, numRadialRegions, clusterFlag, flagPixelwise)
% from showfits.m 11/14/04

% first, read in washins in r,theta from all slices, and then map to x,y positions to gie polar map

% read in all study3.slice*txt files, assume 1 is basal 



% then 

if ~exist('clusterFlag')
   clusterFlag=0;
end
if ~exist('flagPixelwise')
   flagPixelwise=0;
end

[nRegs, nTimes]=size(curves)
bldcurve=curves(1,:)';

nRegs=nRegs-1;
tisscurves=curves(2:nRegs+1,:);
tisscurve=tisscurves';

% read 1st line
for islice=1:numSlices
if (clusterFlag==1)
   filename=strcat(outpath,'fitparams.study',int2str(studyNum),'.slice',int2str(islice),'.clusters','.txt');
elseif (flagPixelwise==1)
   filename=strcat(outpath,'fitparams.study',int2str(studyNum),'.slice',int2str(islice),'.pixelwise','.txt');
else
   filename=strcat(outpath,'fitparams.study',int2str(studyNum),'.slice',int2str(islice),'.',int2str(numAzimuthalRegions),'.',int2str(numRadialRegions),'.txt');
end
fid=fopen(filename ,'r' );
a =fscanf(fid,'%s',3);
a =fscanf(fid,'%s',1);
delta_t=str2num(a(9:14));

indexkwi=1; indexkwo=1; indexfv=1; indext0=1;
while ~isempty(a)
  a =fscanf(fid,'%s',1);
  switch(a)
     case {'flow'}
    	a =fscanf(fid,'%s',1);
    	a =fscanf(fid,'%s',1);   % argh, must be a cleaner way! 
        aa=str2num(a);
        kwi(indexkwi,islice)=aa; indexkwi=indexkwi+1;
     case {'ve'}
    	a =fscanf(fid,'%s',1);
    	a =fscanf(fid,'%s',1);
        aa=str2num(a);
        kwo(indexkwo)=aa;  indexkwo=indexkwo+1;
     case {'fv'}
    	a =fscanf(fid,'%s',1);
    	a =fscanf(fid,'%s',1);
        aa=str2num(a);
        est_fb(indexfv)=aa;  indexfv=indexfv+1;
     case {'t0'}
    	a =fscanf(fid,'%s',1);
    	a =fscanf(fid,'%s',1);
        aa=str2num(a);
        t_delay(indext0)=aa;  indext0=indext0+1;
     otherwise
  end

end


end % slice loop 

%kwi=kwi/2;
est_fb
t_delay
delta_t

% plot circumferential profilses
figure(1); clf; set(gca,'FontSize',16); set(gcf,'Color',[1 1 1]);
hold on
xlabel('Region Number')
ylabel('Flow value')
%for islice=1:sliceNum
      %plot(kwi(nRegs*(islice-1)+1:nRegs*islice),'linewidth',2);
      plot(kwi,'*-','linewidth',2);
avgkwi=mean(mean(kwi));
      plot(avgkwi*ones(1,nRegs),':','linewidth',1);
%end
   Ylim=get(gca,'YLim');
   Ylim(1)=0.0;
   Xlim=get(gca,'Xlim')
   axis([Xlim Ylim]);
legend('Slice 1','Slice 2', 'Slice 3', 'avg')
drawnow

% may want to throw out an outlier from each slice, then scale them to all match in next 3 highest regions for example.

% choose 3 or 4 of the 8 that are most similar:
%for islice=1:numSlices
%   for ireg=1:nRegs
%      for jreg=1:nRegs-1
%        alldiffs(jreg)=kwi(ireg,islice)-kwi(jreg,islice);
%        if ireg==jreg
%          alldiffs(ireg)=-1;
%        end
%      end   
%   end
%end

% better way - look for 3 (or 4) regions most smoothe - least slope, use maxupslopes reoutine??



% figure(3)
% h1=surf(cos(theta1')*rho,sin(theta1')*rho,0*theta1'*rho);

% make polarmaps
   angtorot=-45;
   theta = linspace(0,2*pi,nRegs+1) + angtorot*pi/180;  % don't use last number
   theta=theta(1:nRegs) + 0.5*theta(2);  % first entry zero, want to index center of each of eisght quandrants, hasd problems around zero, 
% arbtitrarily - go to 100 angles:

nX=100; nY=100;
% assume r=slice number? or relate to pixel size and radius of iamge LV?
 x0=nX*0.5+1; 
 y0=nY*0.5;
maxr=sqrt((nX-x0)^2+(nY-y0)^2)+2  % need to recalc if center changes??
% either interp rectoplot to this size or scale
r_spacing=linspace(1,numSlices,maxr);
theta_spacing=linspace(1,nRegs,101);
   theta = linspace(0,2*pi,100+1) + angtorot*pi/180;  % don't use last number
   theta=theta(1:100) + 0.5*theta(2);  % first entry zero, want to index center of each quandrants, hasd problems around zero, 

rectoplot=kwi;
for ireg=1:nRegs
   rectoUpsampled(ireg,:)=interp1(1:numSlices,rectoplot(ireg,:), r_spacing,'cubic');
end
%   rectoUpsampled=interp2(1:nRegs,1:numSlices,rectoplot, 1:nRegs, r_spacing,'cubic');

theta_spacing=linspace(1,3*nRegs,301);
for irad=1:length(r_spacing)
%   tmpp=interp1(1:nRegs,rectoUpsampled(:,irad), theta_spacing,'cubic');
% make wrap around - put 3 in a row, use middle one

   tmppazimuthal=[rectoUpsampled(:,irad); rectoUpsampled(:,irad); rectoUpsampled(:,irad)];
   tmpp=interp1(1:3*nRegs,tmppazimuthal, theta_spacing,'cubic');
   rectoUpsampled2(:,irad)=tmpp(101:200)';
end

maxrad=x0;
polarimg=zeros(nX,nY);
   for y=1:nY
      for x=1:nX
          r=sqrt((x-x0)^2+(y-y0)^2) 
       if r<maxrad
 	  try theta_fromxy=atan2((y-y0),(x-x0))+pi
          catch disp('divide by dzero in atan2')
 	      theta_fromxy=0;
          end
          [value,  closest_theta_index]=sort(abs(theta-theta_fromxy));
% how about linearly interp. along theta, then use nearest neighbor here

%             polarimg(x,y)=(r-floor(r))*rectoUpsampled2(closest_theta_index(1),ceil(r)+1) + (1-(r-floor(r)))*rectoUpsampled2(closest_theta_index(1),floor(r)+1);
             polarimg(x,y)=rectoUpsampled2(closest_theta_index(1),round(r)+1) ;
       end

      end 
   end

figure(2); clf
imglimits=[min(min(kwi)) max(max(kwi))];  % if want max contrast
imagesc(polarimg, imglimits)
colorbar
axis('square')
keyboard

outfilename=strcat(outpath,'fitparams.study',int2str(studyNum),'.slice',int2str(islice),'.',int2str(numAzimuthalRegions),'.',int2str(numRadialRegions),'.polarmap.',int2str(nX),'x',int2str(nY),'.flt');
  ff=fopen(outfilename,'w');
  fwrite(ff,polarimg','float');  % get back origina orientation
  fclose(ff);


return
