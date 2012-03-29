function mpi_resampleCurves(sliceNum, studyNum, outpath, filenamePart1, newRate, skipbeat, original_heartrate, jitter)
%function mpi_resampleCurves(curves, outfilename, newRate, skipbeat, original_heartrate, jitter)
%function upslopes=mpi_upslope(sliceNum,studyNum, outpath, filenamePart1, numSkip

% need to normalize by pre-contrast level  3/1/02

if nargin==0,
    error('sliceNum argument needed for mpi_upslope.');
end


if nargin==0,
    error('args needed for mpi_creaet_rfit');
end
if ~exist('newRate'),
    	newRate=2;
end
if ~exist('original_heartrate'),
    	original_heartrate=60;
end
delta_t=60/original_heartrate;
if (skipbeat)
	delta_t=2*60/original_heartrate;
end
delta_t=delta_t*newRate;

if ~exist('jitter'),
    	jitter=1;   % can be 1..newRate to change which frames deleted
end


curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.mat');
load(curvefilename);
curves = eval(filenamePart1);   % may want to re-do this so not so tricky

%gd_curves =[ bldGds'; tissGds'];  % so in same format as curves

[nRegions, nFrames]=size(curves)
bldcurve=curves(1,:)';
nRegions=nRegions-1;

tisscurves=curves(2:nRegions+1,:);
tisscurve=tisscurves';

% resample
bldcurve=bldcurve(jitter:newRate:end);
tisscurve=tisscurve(jitter:newRate:end,:);
%tisscurve=tisscurve(jitter:newRate:nFrames,:);
nFrames=length(bldcurve);

gdcurves =[ bldcurve'; tisscurve'];  % so in same format as curves
% write out resampled curves... need delta_t too?  
curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.every',int2str(newRate),'.mat');
save(curvefilename, 'gdcurves')


% write out in ascii             %%
% to be added - use only RFIT for now 11/2/02

%/********* Make output rfit file *******************/
outfilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(newRate),'x.rfit');
%outfilename=strcat(outpath,'gdcurves.study',int2str(studyNum),'.slice',int2str
%(sliceNum),'.',int2str(newRate),'x.rfit');

%  fp=fopen(strcat(outfilename,'.gd'),'w')
  fp=fopen(outfilename,'w')
  fprintf(fp,'#%%NTIMES%% %d',nFrames);
  fprintf(fp,'\n#%%NREGIONS%% %d',nRegions+1);
  fprintf(fp,'\n#%%NCORREL%%\t%d',1);
  fprintf(fp,'\n#%%CORLIST%%\t%d',1);
  fprintf(fp,'\n#%%TRUE_UNCERT%%\t%d',1);

  fprintf(fp,'\n#\tStart\tStop\tActive');

  acq_time=0.1;
%  time_between_frames=1.666;  % assumes 1.2 beats/sec (hr=72, 1 beat=833msec)
  time_between_frames=delta_t;  % 2.3;  % 12/4/01 for DG P062801 assumes 50 beats/min
  for i=0:nFrames-1
     fprintf(fp,'\n%f\t%f\t%f',i*time_between_frames,i*time_between_frames+acq_time,acq_time);
  end
  fprintf(fp,'\n\n');

  fprintf(fp,'# Region 1');
  fprintf(fp,'\n#%%LABEL%%\tBlood Input');
  fprintf(fp,'\n#cts/vox/acq_time\tUncert\tInput Correlations');

  tim=acq_time;
  for iFrame=1:nFrames
       input(iFrame) = bldcurve(iFrame);
       if (input(iFrame) < 0) 
	   input(iFrame)=0;
	end
%       std=sqrt(input(iFrame));  	% Poisson noise
       std=0; 			 	% arbitrary Gaussian noise
       if (std == 0.0) 
	   std = 1.0 ;
       end
       fprintf(fp,'\n%f\t%f\t%f',input(iFrame),std,1.0);
       tim=tim+acq_time;
  end

  for ireg=1:nRegions
     fprintf(fp,'\n\n');
     fprintf(fp,'# Region %d',ireg+1);
     fprintf(fp,'\n#%%LABEL%%\tTissue');
     fprintf(fp,'\n#cnts/vox/acq_time\tUncert');

     tim=acq_time;
     for iFrame=1:nFrames
        uptake(iFrame) = tisscurve(iFrame,ireg) ;
       if (uptake(iFrame) < 0) 
	   uptake(iFrame)=0;
	end
%        std=sqrt(uptake(iFrame));	% Poisson noise
        std=10; 		 	% arbitrary Gaussian noise
        if (std == 0.0) 
            std = 1.0;
        end
        fprintf(fp,'\n%f\t%f\t%f',uptake(iFrame),std,0.0);
        tim=tim+acq_time;
     end
   end


% /* Add some carriage returns for RFIT */
fprintf(fp,'\n\n\n');

fclose(fp)
