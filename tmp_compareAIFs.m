
outpath='Output/'
seriesNum=12
sliceNum=1
curvefilename=strcat(outpath,'curves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat')
load(curvefilename)
bld1=curves(1,:);
% showcurves(curves,'Sig. Intensity')


sliceNum=2
curvefilename=strcat(outpath,'curves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat')
load(curvefilename)
bld2=curves(1,:);


sliceNum=3
curvefilename=strcat(outpath,'curves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat')
load(curvefilename)
bld3=curves(1,:);


clf; hold on;
plot(bld1,'b')
plot(bld2,'r')
plot(bld3,'g')
legend('bld1', 'bld2', 'bld3')
