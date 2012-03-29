function curves=mpi_useAlternateAIF(sliceNum,studyNum, outpath, numPreContrast_bld, numPreContrast_tiss, numSkip, scale_to_AIF)

numCCs=12;
% if pre-bolus has fewer time frames just replace time frames up to that point...  7/5/05

% note: doesn't care if pixelwise or not, just read in a file of curves

disp('WARNING!!!  THIS WILL OVERWRITE THE CURVES FILES TO PUT IN ALTERNATE AIF. SCALING MAY BE OFF!!!  Edit mpi2d if needed.')
disp('Hit CNTRL-C to abort if you do not want curves overwritten! Easy to get back by re-running 3.2 ')

studyNumAIF=input('Enter study or series number to get AIF from \n')
sliceNumAIF=input('Enter slice number to get AIF from \n')

% not logging this, figure can check which matches which if necessary by doing diffs and looking at dates...

% read in alternate AIF 
   curvefilename=strcat(outpath,'curves.study',int2str(studyNumAIF),'.slice',int2str(sliceNumAIF),'.mat');
   try load(curvefilename);
   catch disp('cant seem to find curves.study file to load, may not have run 3.2 yet')
   end
   aif_curve=curves(1,:)';
   [nRegsAIF, nTimesAIF]=size(curves);

   clear curves

   curvefilename=strcat(outpath,'curves.study',int2str(studyNum),'.slice',int2str(sliceNum),'.mat');
   load(curvefilename);

   si_curves = curves;
   [nRegs, nTimes]=size(curves);
   nRegs=nRegs-1;
   bldcurve=si_curves(1,:)';

   % scale pre-contrast to match in blood:
   init_aif=sum(aif_curve(numSkip+1:numPreContrast_bld+numSkip))/numPreContrast_bld;
   init_bld=sum(bldcurve(numSkip+1:numPreContrast_bld+numSkip))/numPreContrast_bld;
%   bldcurve=aif_curve*10*init_bld/init_aif;

 
   disp('USING ALTERNATE AIF')
if (scale_to_AIF==1) 
% don't scale up by 10, just scale resulting param. ktrans by 10...
   bldcurve(1:nTimesAIF)=aif_curve*init_bld/init_aif;
   disp('scaling aif by init_bld/init_aif NOT times 10')
   init_bld/init_aif
else
   bldcurve(1:nTimesAIF)=aif_curve*numCCs;
%   bldcurve(1:nTimesAIF)=aif_curve*4;
end

   tisscurves = si_curves(2:nRegs+1,:);
   tisscurve=tisscurves';



curves =[ bldcurve'; tisscurve'];  % so in same format as curves


showcurves([bldcurve'; tisscurve'],'Alternate AIF shown');
return;
