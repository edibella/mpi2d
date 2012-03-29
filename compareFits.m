function compareFits(file1, file2, sliceNum, studyNum,outpath)
%COMPAREFITS (file1, file2, sliceNum, studyNum,outpath)
% note hard-coded nrows, ncols
% don't need last 3 args if filesizes same)
% does some cropping, and percent differences and plots it

nrows=144; ncols=92;


if nargin==0,
    error('arguments needed for mpi_remakeOriginalWithFits');
 end

% read in data
%filename=strcat('fitparams.study',int2str(studyNum),'.slice',int2str(sliceNum),'.txt');
fid=fopen(file1 ,'r' );
img1=fread(fid,[nrows,ncols],'float');
fclose(fid)
img1=img1';
imagesc(img1)

% read in data
fid=fopen(file2 ,'r' );
img2=fread(fid,[nrows,ncols],'float');
img2=img2';
imagesc(img2)


tmpImg=zeros(nrows,ncols);
[X,Y,xCenter, yCenter, start_angle,bw1, bw3]= mpi_loadContours(tmpImg,sliceNum,studyNum, outpath);

% find max and min X and Y to set crop boundaries:
   maxX=max(X)+1; minX=min(X)-1;
   maxY=max(Y)+1; minY=min(Y)-1;
%   imgParams=imgParams(minY:maxY,minX:maxX,:);

% embed (zero-pad) the original:
img1new=zeros(nrows,ncols);
img1new(minY:maxY,minX:maxX)=img1(:,:);
img2new=zeros(nrows,ncols);
img2new(minY:maxY,minX:maxX)=img2(:,:);


img1=img1new;
img2=img2new;
imagesc(img2)
pause


nX=length(X);

for j = 1 : nX
    kwi1(j)=img1(Y(j), X(j));
    kwi2(j)=img2(Y(j), X(j));
    diff_kwi(j)=img1(Y(j), X(j))-img2(Y(j), X(j));
end


%end

%diff_kwi=kwi1 - kwi2(index); 

clf
hold on
set(gca,'FontSize',16)
%plot(diff_kwi,'xm','Linewidth',5)
plot(100*diff_kwi./kwi1,'xm','Linewidth',5)
title('Differences in washin parameter')
xlabel('Pixel number')
ylabel('Percent Difference')

length(diff_kwi)
meandiff=mean(abs(diff_kwi))
maxdiff=max(abs(diff_kwi))
sd_diff=std(abs(diff_kwi))


allkwi=[kwi1 kwi2];
meanallkwi=mean(allkwi)
disp('Percent differences:')
meanpercentdiff=mean(abs(diff_kwi))/mean(allkwi)
maxpercentdiff=max(abs(diff_kwi))/mean(allkwi)
sd_percentdiff=std(abs(diff_kwi))/mean(allkwi)

meanpercentdiff=100*mean(abs(diff_kwi)./kwi1)
maxpercentdiff=100*max(abs(diff_kwi)./kwi1)
sd_diff=std(abs(diff_kwi))
sd_percentdiff=100*std(abs(diff_kwi)./kwi1)


figure(2)
plot( kwi1, kwi2,'x','Linewidth',3)
hold on
% regreession:

coefs=polyfit(kwi1,kwi2,1)   % first-order polynomial (line int. and slope)
slope=coefs(1)
newy=polyval(coefs,kwi1);
%lsq_dist=sum(line-data)

plot( kwi1, newy,'k','Linewidth',2)
set(gca,'FontSize',16)
set(gcf,'Color',[1 1 1])   % to get rid of grey border
xlabel('Washin truth')
ylabel('Washin estimate')



keyboard
% why did this stuff stop working? non-invertible with inv()
%twoplot(kwi2, kwi1)

%[tmpcoeffs,Serr] = lregress(kwi2, kwi1)
%          slope=tmpcoeffs(2,1);




return;




%function index=computeIndex(region,tmpImg,sliceNum,studyNum,outpath)

%return;

