function mpi_viewResultsEndoAndEpi(studyNum, sliceVector, outpath, numAzimuthalRegions, numRadialRegions, clusterFlag, flagPixelwise)


% HARD CODED FOR CVTMRF    12/6/04



if ~exist('clusterFlag')
   clusterFlag=0;
end
if ~exist('flagPixelwise')
   flagPixelwise=0;
end

nRegs=numAzimuthalRegions*numRadialRegions;


% customize to make apex to base:     (tjhis is for cvtmrf)
indexslice=0;
numSlices=length(sliceVector);
for islice=sliceVector;   %[2 3 1]

if (clusterFlag==1)
   filename=strcat(outpath,'fitparams.study',int2str(studyNum),'.slice',int2str(islice),'.clusters','.txt');
elseif (flagPixelwise==1)
   filename=strcat(outpath,'fitparams.study',int2str(studyNum),'.slice',int2str(islice),'.pixelwise','.txt');
else
   filename=strcat(outpath,'fitparams.study',int2str(studyNum),'.slice',int2str(islice),'.',int2str(numAzimuthalRegions),'.',int2str(numRadialRegions),'.txt')
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

%kwi=kwi/2;
%Hct=0.45;
%kwo=(1-Hct)*kwi./kwo;
est_fb
t_delay
delta_t

kwiEndo=kwi(1:numAzimuthalRegions,:);
kwiEpi=kwi(numAzimuthalRegions+1:2*numAzimuthalRegions,:);
jj=1;
for ii=1:2:length(kwiEndo)
   kwi(ii)=kwiEndo(jj);
   kwi(ii+1)=kwiEpi(jj+1);
   jj=jj+1;
end

% plot circumferential profilses
studyDescriptor='CVTMRN ';
vectorDescriptor='EndoEpiInterleaved ';
figure(2)
if studyNum==10
  colorVector=['b','g','r','c','m']
  lineType='--'
else
  colorVector=['b','g','r','c','m']
%  colorVector=['m','y','k','b','g']
  lineType='-*'
end
plotCprof(kwi, studyNum, vectorDescriptor, studyDescriptor, lineType, colorVector)  % last parameter is optional


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



% h1=surf(cos(theta1')*rho,sin(theta1')*rho,0*theta1'*rho);

% make polarmaps
%makePolar(outpath, studyNum, numAzimuthalRegions, numRadialRegions/2, numSlices, kwiEndo )

return
