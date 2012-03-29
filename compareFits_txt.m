%function compareFits(file1, file2)
%function compareFits_txt(file1, file2, sliceNum, studyNum,outpath)
function [params timecurves]=compareFits_txt(sliceNum, studyNum, sliceNum2, studyNum2)

outpath='Output/'

% don't need last 3 args if filesizes same)


if nargin==0,
    error('arguments needed for compareFits_txt');
 end
 

% read in data
%filename=strcat('fitparams.study',int2str(studyNum),'.slice',int2str(sliceNum),'.txt');

numAzimuthalRegions=6;
numRadialRegions=1;
%numAzimuthalRegions=16;

studyNumAIF=studyNum+99
sliceNumAIF=sliceNum
scaleAIF=10


filename=strcat(outpath,'flowvalues.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions),'.',int2str(numRadialRegions),'.txt') 
filename=strcat(outpath,'flowvalues.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions),'.',int2str(numRadialRegions),'_fixedDelay0.AIF_',int2str(studyNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.txt.full.txt') 

fid=fopen(filename ,'r' );
a =fscanf(fid,'%s',3);
a =fscanf(fid,'%s',1);
delta_t=str2num(a(9:14));

indexkwi=1; indexkwo=1; indexfv=1; indext0=1;
while ~isempty(a)
  a =fscanf(fid,'%s',1);
  switch(a)
     case {'Ktrans'}
        a =fscanf(fid,'%s',1);
        a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
        aa=str2num(a);
        kwi(indexkwi)=aa; indexkwi=indexkwi+1;
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

flow1=kwi;
kwi=kwi/2;
Hct=0.45;
kwo=kwi./((1-Hct)*kwo);


nSect=length(kwi);
fclose(fid);

%   iseries=first_series+(islice-1)*1000+isub
% studyNum=first_series+(sliceNum-1)*1000+isub

curvefilename=strcat(outpath,'deltaSIcurves.study',int2str(studyNum),'.slice',int2str(sliceNum),'.AIF_',int2str(studyNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.mat') 

load(curvefilename)
bld1=deltaSIcurves(1,:);
deltaSIcurves1=deltaSIcurves;
clear deltaSIcurves;

%check curves before, should be no time stamps and no special AIF
curvefilename=strcat(outpath,'curves.study',int2str(studyNum),'.slice',int2str(sliceNum),'.mat')
load(curvefilename)
curves1=curves;
clear curves;


studyNumAIF=studyNum2-2
sliceNumAIF=sliceNum2
%scaleAIF=10

filename=strcat(outpath,'fitparams.study',int2str(studyNum2),'.slice',int2str(sliceNum2),'.',int2str(numAzimuthalRegions),'.',int2str(numRadialRegions),'.txt') 
filename=strcat(outpath,'flowvalues.study',int2str(studyNum2),'.slice',int2str(sliceNum2),'.',int2str(numAzimuthalRegions),'.',int2str(numRadialRegions),'_fixedDelay0.AIF_',int2str(studyNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.txt.full.txt') 

fid=fopen(filename ,'r' );
a =fscanf(fid,'%s',3);
a =fscanf(fid,'%s',1);
delta_t2=str2num(a(9:14));

indexkwi=1; indexkwo=1; indexfv=1; indext0=1;
while ~isempty(a)
  a =fscanf(fid,'%s',1);
  switch(a)
     case {'Ktrans'}
        a =fscanf(fid,'%s',1);
        a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
        aa=str2num(a);
        kwi2(indexkwi)=aa; indexkwi=indexkwi+1;
     case {'ve'}
        a =fscanf(fid,'%s',1);
        a =fscanf(fid,'%s',1);
        aa=str2num(a);
        kwo2(indexkwo)=aa;  indexkwo=indexkwo+1;
     case {'fv'}
        a =fscanf(fid,'%s',1);
        a =fscanf(fid,'%s',1);
        aa=str2num(a);
        est_fb2(indexfv)=aa;  indexfv=indexfv+1;
     case {'t0'}
        a =fscanf(fid,'%s',1);
        a =fscanf(fid,'%s',1);
        aa=str2num(a);
        t_delay2(indext0)=aa;  indext0=indext0+1;
     otherwise
  end

end

flow2=kwi2;
kwi2=kwi2/2;
Hct=0.45;
kwo2=kwi2./((1-Hct)*kwo2);

%Output/deltaSIcurves.study54092.slice3.AIF_54090_3_10.mat
curvefilename=strcat(outpath,'deltaSIcurves.study',int2str(studyNum2),'.slice',int2str(sliceNum2),'.AIF_',int2str(studyNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.mat') 

%curvefilename=strcat(outpath,'curves.study',int2str(studyNum2),'.slice',int2str(sliceNum2),'.mat')
load(curvefilename)
bld2=deltaSIcurves(1,:);
deltaSIcurves2=deltaSIcurves;
clear deltaSIcurves;


%check curves before, should be no time stamps and no special AIF
curvefilename=strcat(outpath,'curves.study',int2str(studyNum2),'.slice',int2str(sliceNum2),'.mat')
load(curvefilename)
curves2=curves;
clear curves;



params.flow1=flow1;
params.flow2=flow2;
timecurves.study1=deltaSIcurves1;
timecurves.study2=deltaSIcurves2;

% global scale:
%flow1=flow1./sum(flow1)
%flow2=flow2./sum(flow2)

diff_flow=flow1 - flow2

clf
hold on
set(gca,'FontSize',16)
%plot(diff_kwi,'xm','Linewidth',5)
plot(100*diff_flow./flow1,'xm','Linewidth',5)
title('Differences in flow est. parameter')
xlabel('Pixel number')
ylabel('Percent Difference')

pause


diff_kwi=diff_flow;
length(diff_kwi)
meandiff=mean(abs(diff_kwi))
maxdiff=max(abs(diff_kwi))
sd_diff=std(abs(diff_kwi))


kwi1=kwi;
clf
kwi1=flow1;
kwi2=flow2;

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
xlabel('Washin first file')
ylabel('Washin 2nd filename')





figure; clf; hold on
plot(bld1,'r')
plot(bld2);
plot(deltaSIcurves1(2:7,:)','r')
plot(deltaSIcurves2(2:7,:)','b')
legend('first AIF ', '2nd AIF')

% figure; clf; hold on
% plot(curves1','m')
% plot(curves2','y')
	 
mean(params.flow1)
std(params.flow1)
mean(params.flow2)
std(params.flow2)

% why did this stuff stop working? non-invertible with inv()
%twoplot(kwi2, kwi1)

%[tmpcoeffs,Serr] = lregress(kwi2, kwi1)
%          slope=tmpcoeffs(2,1);




return;




%function index=computeIndex(region,tmpImg,sliceNum,studyNum,outpath)

%return;

