function mpi_remakeClusters(img,sliceNum, studyNum, outpath )


% idea is to read in cluster image, embed in standard size, mask, 
%then read in fits, and write out cluster image with washin values. Size will be smaller than standard 93x79, this is a problem... need to reverse crops...
% NO, just apply crops to other estimates

if nargin==0,
    error('arguments needed for mpi_remakeClusters');
 end


[X,Y,xCenter, yCenter, start_angle,bw1, bw3]= mpi_loadContours(img,sliceNum, studyNum,outpath);
% find max and min X and Y to set crop boundaries:
   maxX=max(X)+1; minX=min(X)-1;
   maxY=max(Y)+1; minY=min(Y)-1;


% read in data
nClusters=9;
nFrames=50;
nrows=maxX-minX+1; ncols=maxY-minY+1;
filename=strcat(outpath,'croppedForClustering.',int2str(maxX-minX+1),'x',int2str(maxY-minY+1),'x',int2str(nFrames),'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(nClusters),'.tacs.img')
%[PATH,NAME,EXT,VER]= fileparts(filename);
%data=load(filename);

fid=fopen(filename,'r');
img=fread(fid, [nrows ncols],'float');
clf
figure(1)

imagesc(img');
axis('square')


% read 1st line
filename=strcat('fitparams.study',int2str(studyNum),'.slice',int2str(sliceNum),'.txt');
fid=fopen(filename ,'r' );
a =fscanf(fid,'%s',3);
a =fscanf(fid,'%s',1);
delta_t=str2num(a(9:14));

indexkwi=1; indexkwo=1; indexfv=1; indext0=1;
while ~isempty(a)
  a =fscanf(fid,'%s',1);
  switch(a)
     case {'kwi'}
      a =fscanf(fid,'%s',1);
        aa=str2num(a);
        kwi(indexkwi)=aa; indexkwi=indexkwi+1;
     case {'kwo'}
      a =fscanf(fid,'%s',1);
        aa=str2num(a);
        kwo(indexkwo)=aa;  indexkwo=indexkwo+1;
     case {'fv'}
      a =fscanf(fid,'%s',1);
        aa=str2num(a);
        est_fb(indexfv)=aa;  indexfv=indexfv+1;
     case {'t0'}
      a =fscanf(fid,'%s',1);
        aa=str2num(a);
        t_delay(indext0)=aa;  indext0=indext0+1;
     otherwise
  end

end


nClusters=length(kwi);
nSect=nClusters;
nX     = length(X);
imgParams=zeros(nrows,ncols);


for i=1:nClusters
   clear x y
   [x, y]=find(img==(i-1));
   for j=1:length(x)
      imgParams(x(j),y(j))=kwi(i);
   end
end

imagesc(imgParams')

%keyboard


outfilename=strcat(outpath,'img.study',int2str(studyNum),'.slice',int2str(sliceNum),'newOriginal.clusters.washins.float');

ff=fopen(outfilename,'w')
fwrite(ff,imgParams,'float');
fclose(ff);

return;
