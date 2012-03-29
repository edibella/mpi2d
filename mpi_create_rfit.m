function mpi_create_rfit(curves, outfilename, delta_t)

if nargin==0,
    error('args needed for mpi_creaet_rfit');
end
if ~exist('delta_t'),
   delta_t=2
end

[nRegions, nFrames]=size(curves)
bldcurve=curves(1,:)';
nRegions=nRegions-1;

tisscurves=curves(2:nRegions+1,:);
tisscurve=tisscurves';


%/********* Make output rfit file *******************/
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
