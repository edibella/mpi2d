
outpath='Output/'
first_series=52000  %25000  %53000

figure; clf; hold on;

for islice=3:5  %2:8
    
    isub=99
    isub=92
   % isub=0
   iseries=first_series+(islice-1)*1000+isub

curvefilename=strcat(outpath,'curves.study',int2str(iseries),'.slice',int2str(islice),'.mat')
load(curvefilename)
bld=curves(1,:);

plot(bld);
drawnow
pause
end


%bld(islice,:)=curves(1,:);
% showcurves(curves,'Sig. Intensity')

% 
% sliceNum=2
% curvefilename=strcat(outpath,'curves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat')
% load(curvefilename)
% bld2=curves(1,:);
% 
% 
% sliceNum=3
% curvefilename=strcat(outpath,'curves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat')
% load(curvefilename)
% bld3=curves(1,:);




% clf; hold on;
% plot(bld1,'b')
% plot(bld2,'r')
% plot(bld3,'g')
legend('bld1', 'bld2', 'bld3')
