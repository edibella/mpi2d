function precontrastT1=fit_subsets(ProcessingDir, seriesnum, slicenum, seriesNumOffset,precontrastT1)  
%precontrastT1 got moved to last param to allow for ganeshes' modified series number soffset 
% last parameter optional. Will fit 3 curves, write out all results, even mpi2d flows
%  as in fit_subsets('P040110',18000,1)
% from /v/raid1/npack/Code/ReconCode/BatchFile_AIFsubsets_123.m

raid1path=[ProcessingDir '/Output/'];
if(~exist('seriesNumOffset'))
    seriesNumOffset = 0;
end
series=seriesnum;
slice=slicenum; 
%series=18000;slice=1;

 for isubset=(1:3)+seriesNumOffset
   sfile=strcat(raid1path,'curves.study',int2str(series+isubset),'.slice',int2str(slice),'.mat');      % 1st subset corresponds to lowest SRT (ie. 1st signal)
   load (sfile)
   uload(:,isubset)=curves(1,:);
 end
 
 % now swap them so 2nd subset is highest SRT!  % fixed calc_T1 5/20/10 so
 % don't need to do this   EVRD
%  tmpp=uload(:,2);
%  uload(:,2)=uload(:,3);
%  uload(:,3)=tmpp;

threeCurves=uload;  

% first compute T1_0 for first injection:
%precontrastT1=

if exist('precontrastT1')==1
   [Gdconc Params]=calc_T1_subset_123_72ray(threeCurves,slice,precontrastT1)
else
    [Gdconc Params]=calc_T1_subset_123_72ray(threeCurves,slice)
    precontrastT1=Params.EstPrecontrast
end

outfile=strcat(raid1path,'MultiSRT_Vals_M0_fA_T1_chi2err.series',int2str(series),'.slice',int2str(slice),'.mat');
save (outfile,'Params')
%%% The lines of code below were adapted from the original file "gdAIF.m" in order to 
%%% compute the re-write the [Gd] curves  into the output directory for each patient.
loadfilename=strcat(raid1path,'deltaSIcurves.study',int2str(series),'.slice',int2str(slice),'.mat');
savefilename=strcat(raid1path,'deltaSIcurves.study',int2str(series),'.slice',int2str(slice+10),'.mat');
load (loadfilename); deltaSIcurves(1,:)=gdplot(Gdconc,deltaSIcurves); save (savefilename,'deltaSIcurves')

cd(ProcessingDir);
% !cp series23.slice3.par mpi2d.par 
% infile=strcat('series',int2str(seriesnum),'.slice',int2str(slice),'.par');
% copyfile(infile,'mpi2d.par')

%addpath('/v/raid1/npack/Code/mpi2dCode_THK_092209/')
%makeParFileED(infolder, infile, outfolder, templateParFile, slice, subset, study,tmax)
%    subset=1
% seriesnum=24*1000+subset
% infolder=strcat('/v/raid1/ed/MRIdata/Cardiac/Trio/P030210/ReconData/CV_Radial7Off_triple_3.1ml(',int2str(seriesnum),')/')
% parfile= makeParFile_forSubsets(infolder, outfolder, templateParFile, slice, subset, seriesnum, copyFromSeriesNum)

%addpath('/Users/ed/Documents/ShareWithWindows/matlab/generalCodes')
makeParFileForGadAIF(series, slice, slice+10)   %puts into mpi2d.par 
%parfile=strcat('rad.series',int2str(seriesnum),'.slice',int2str(slice),'.subset',int2str(subset),'.par')
%copyfile(parfile, 'mpi2d.par')

mpi2d scaleFrame1=40 scaleFrame2=42 stg=[4.1]  %[4.2032]




return







% disp(' will switch directories to compare to TCR results ')
% pause
% 
% 
%  raid1path='/v/raid1/npack/Processing_72ray_subset123/';
% % Patient='P070909';
%  
%    sfile=strcat(raid1path,Patient,'_nufft_24rays/Output_subset1/curves.study',int2str(series),'.slice',int2str(slice),'.mat');      % 1st subset corresponds to lowest SRT (ie. 1st signal)
%    load (sfile)
%    uload(:,1)=curves(1,:);
%    sfile=strcat(raid1path,Patient,'_nufft_24rays/Output_subset2/curves.study',int2str(series),'.slice',int2str(slice),'.mat');      % 2nd subset corresponds to highest SRT (ie. 3rd signal)
%    load (sfile)
%    uload(:,2)=curves(1,:);
%    sfile=strcat(raid1path,Patient,'_nufft_24rays/Output_subset3/curves.study',int2str(series),'.slice',int2str(slice),'.mat');      % 3rd subset corresponds to medium SRT (ie. 2nd signal)
%    load (sfile)
%    uload(:,3)=curves(1,:);
%    
% 
% 
% 
% %cd /v/raid1/npack/Processing_72ray_subset123/P081309_nufft_24rays/
% % from /v/raid1/npack/Code/ReconCode/BatchFile_AIFsubsets123.m
% series=21;slice=3;precontrastT1=1664  %767;
% calc_T1_subset_123_eq(Patient,uload,slice,precontrastT1)
% !cp series23.slice3.par mpi2d.par                   
% addpath('/v/raid1/npack/Code/mpi2dCode_THK_092209/')
% mpi2d sliceNumAIF=13 scaleFrame1=40 scaleFrame2=42 stg=[4.2032]

%calc_T1_subset_123_eq(Patient,threeCurves,slice,precontrastT1)
