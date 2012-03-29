function cine_curv_fit(vart,noi)
load Y;
load X;
load init_bld;
load init_tiss;
load tiss_avgSpatial;

if(vart==1)
load varg; % fmincon -non-linear
end

if(vart==0)
load vare; % murase
end

if(vart==1)
est_cur=varg; 
else
    est_cur=vare;
end

disp('Done loading of variables')


nRegs=size(est_cur,2);

for i=1:nRegs
 orig_tisscurve(:,i)=(est_cur(:,i)*init_tiss(i))/(tiss_avgSpatial);
 orig_tisscurve(:,i)=orig_tisscurve(:,i)+init_tiss(i);
end
 outpath='Output/';
% 
infilename=strcat(outpath,'cinemri1.study',int2str(studyNum),'.slice',int2str(sliceNum),'.mat');
load(infilename);

cinemri=cinemri1;


img=cinemri;
new_tisscurves=orig_tisscurve';

cinemri_curv_fit=img;
nX=length(X);
not_fit=200*ones(1,size(new_tisscurves,2));
counter=0;

for j = 1 : nX
    a=find(new_tisscurves(j,:)<10000);
    if(a)
      cinemri_curv_fit(Y(j), X(j), :)= new_tisscurves(j,:);
  else
      counter=counter+1;
      cinemri_curv_fit(Y(j), X(j), :)=not_fit;
    end
end
disp('Total no.of pixels')
nRegs

disp('No.of pixels not fit')
counter     

outfilename=strcat(outpath,'cinemri_curv_fit.study',int2str(noi),'_',int2str(studyNum),'.slice',int2str(sliceNum),'.mat');
save(outfilename, 'cinemri_curv_fit')
