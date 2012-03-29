function mpi_viewResults(studyNum, sliceVector, outpath, numAzimuthalRegions, numRadialRegions, clusterFlag, flagPixelwise)
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

nRegs=numAzimuthalRegions*numRadialRegions+1;

nRegs=nRegs-1;

% customize to make apex to base:     (tjhis is for cvtmrf)
indexslice=0;
numSlices=length(sliceVector);
for islice=sliceVector
%for islice=[3]

%for islice=1:numSlices
if (clusterFlag==1)
   filename=strcat(outpath,'fitparams.study',int2str(studyNum),'.slice',int2str(islice),'.clusters','.txt')
elseif (flagPixelwise==1)
   filename=strcat(outpath,'fitparams.study',int2str(studyNum),'.slice',int2str(islice),'.pixelwise','.txt')
else
   filename=strcat(outpath,'fitparams.study',int2str(studyNum),'.slice',int2str(islice),'.',int2str(numAzimuthalRegions),'.',int2str(numRadialRegions),'.txt')
%   filename=strcat(outpath,'flowvalues.study',int2str(studyNum),'.slice',int2str(islice),'.',int2str(numAzimuthalRegions),'.',int2str(numRadialRegions),'.AIF_9_1_10.txt')
%   filename=strcat(outpath,'flowvalues.study',int2str(studyNum),'.slice',int2str(islice),'.',int2str(numAzimuthalRegions),'.',int2str(numRadialRegions),'.AIF_20_1_14.txt')
end
fid=fopen(filename ,'r' );
a =fscanf(fid,'%s',3);
a =fscanf(fid,'%s',1);
delta_t=str2num(a(9:14));

indexkwi=1; indexkwo=1; indexfv=1; indext0=1;
indexslice=indexslice+1;
while ~isempty(a)
  a =fscanf(fid,'%s',1);
  switch(a)
     case {'flow'}
    	a =fscanf(fid,'%s',1);
    	a =fscanf(fid,'%s',1);   % argh, must be a cleaner way! 
        aa=str2num(a);
        kwi(indexkwi,indexslice)=aa; indexkwi=indexkwi+1;
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

%kwi=kwi/2;   % not considering plasma, but whole blood flow?
est_fb
t_delay
delta_t

% plot circumferential profilses
figure(studyNum); clf; set(gca,'FontSize',16); set(gcf,'Color',[1 1 1]);
hold on
xlabel('Region Number')
ylabel('Flow value')
      plot(kwi,'*-','linewidth',2);
avgkwi=mean(mean(kwi));
% maybe better, compute avg of highest half:
for islice=0:numSlices-1
  for jj=1:nRegs 
  aa(jj+islice*nRegs)=kwi(jj,islice+1);
  end
end
  bb=sort(aa);
  avgkwi=mean(bb(round(0.5*length(bb)):length(bb)));
      plot(avgkwi*ones(1,nRegs),':','linewidth',1);
%end
   Ylim=get(gca,'YLim');
   Ylim(1)=0.0;
   Xlim=get(gca,'Xlim')
   axis([Xlim Ylim]);
legend('Slice 1 - most apical','Slice 2', 'Slice 3 - most basal', 'avg of higher half')
%title(strcat('CVTMRF Series ',int2str(studyNum)))
%title(strcat(studyTitle,' Series ',int2str(studyNum)))
title(strcat(' Series ',int2str(studyNum)))
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
%makePolar(outpath, studyNum, numAzimuthalRegions, numRadialRegions, numSlices, kwi);

return
