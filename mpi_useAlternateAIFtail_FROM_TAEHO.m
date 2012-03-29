function deltaSIcurves=mpi_useAlternateAIFtail(seriesNumAIF, sliceNumAIF, scaleAIF, sliceNum,seriesNum, outpath, num_bld1a, num_bld2a, numSkip)

numCCs=scaleAIF;
% if pre-bolus has fewer time frames just replace time frames up to that point...  7/5/05


%disp('WARNING!!!  THIS WILL OVERWRITE THE CURVES FILES TO PUT IN ALTERNATE AIF. SCALING MAY BE OFF!!!  Edit mpi2d if needed.')
%disp('Hit CNTRL-C to abort if you do not want curves overwritten! Easy to get back by re-running 3.2 ')
%studyNumAIF=input('Enter study or series number to get AIF from \n')
%sliceNumAIF=input('Enter slice number to get AIF from \n')

% read in alternate AIF 
   curvefilename=strcat(outpath,'deltaSIcurves.study',int2str(seriesNumAIF),'.slice',int2str(sliceNumAIF),'.mat')
   try load(curvefilename);
   catch disp('cant seem to find curves.series file to load, may not have run 3.2 yet')
   end
   aif_curve=deltaSIcurves(1,:)';
   [nRegsAIF, nTimesAIF]=size(deltaSIcurves);

 first_deltaSI=deltaSIcurves;
   clear deltaSIcurves

   curvefilename=strcat(outpath,'deltaSIcurves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat');
   load(curvefilename);
   second_deltaSI=deltaSIcurves;

   si_curves = deltaSIcurves;
   [nRegs, nTimes]=size(deltaSIcurves);
   nRegs=nRegs-1;
   bldcurve=si_curves(1,:)';
   
  
  
   % scale last 5 points of AIF to match?
   init_aif=mean(aif_curve(num_bld1a:num_bld2a));%thk
   init_bld=mean(bldcurve(num_bld1a:num_bld2a));%thk
%   bldcurve=aif_curve*10*init_bld/init_aif;

% % if   sliceNumAIF >10 % by thk
% %     tailscale=1;
% % else
% %     tailscale=init_bld/init_aif; 
% % end
 tailscale=init_bld/init_aif; 

%%% To scale the AIF from this slice to another slice(since above the mean
%%% pre-contrast frames are near 0, I'll use the original curves files to get the pre-contrast scale factor....
%%% NP 04/27/07 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tmp1=strcat(outpath,'curves.study',int2str(seriesNumAIF),'.slice',int2str(sliceNumAIF),'.mat');
% load(tmp1);
% tmpAIFcurve=curves(1,:)';
% tmp2=strcat(outpath,'curves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat');
% load(tmp2);
% tmpBLDcurve=curves(1,:)';
% init_aif=sum(tmpAIFcurve(numSkip+1:numPreContrast_bld+numSkip))/numPreContrast_bld;
% init_bld=sum(tmpBLDcurve(numSkip+1:numPreContrast_bld+numSkip))/numPreContrast_bld;
% %%% to do this, scale_to_AIF should be set to 1 below!
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   disp('USING ALTERNATE AIF')
%scale_to_AIF=1;  % this is set to 1 only if scaling AIF from another slice
scale_to_AIF=1;  % this is set to 1 only if scaling AIF from another slice
add_inputsFlag=0;
summedAIF=0;

if (scale_to_AIF==1) 
% don't scale up by 10, just scale resulting param. ktrans by 10...
  % bldcurve(1:nTimesAIF)=aif_curve*init_bld/init_aif;
   bldcurve(1:nTimesAIF)=aif_curve*abs(tailscale); %not negative by thk
   disp('scaling aif by init_bld/init_aif NOT times 10')
   %init_bld/init_aif
   abs(init_bld/init_aif) %thk
else
   bldcurve(1:nTimesAIF)=aif_curve*numCCs;
   if add_inputsFlag==1
      delta_t=1.0463;
      disp('USING HARD CODED DELTA_T HERE!!!!')
      for k=1:numCCs-1
        t_delay=-k/6;   % cc/sec
        t_delay = t_delay/delta_t;   % convert from seconds to frames
        tmpbldcurve(1-fix(t_delay):nTimes)=bldcurve(1:nTimes+fix(t_delay));
        tmpbldcurve(1:1-fix(t_delay))=0;
        t_delay=t_delay-fix(t_delay);
        tmpbldcurve = interp1(1:nTimes,tmpbldcurve,1+t_delay:nTimes+t_delay,'splines',0);
        summedAIF=tmpbldcurve + summedAIF;
      end
      bldcurve(1:nTimesAIF)=summedAIF(1:nTimesAIF);
    end

end
   tisscurves = si_curves(2:nRegs+1,:);
   tisscurve=tisscurves';

   %%%% In both bldcurve lines below, only take range from 1:nTimes of low-
   %%%% dose AIF so lengths are the same.  This won't work if the length of
   %%%% the low-dose AIF is LESS than the new high-dose AIF   NP 10/27/06
deltaSIcurves =[ bldcurve(1:nTimes)'; tisscurve'];  % so in same format as curves

%keyboard
showcurves([bldcurve(1:nTimes)'; tisscurve'],'Alternate AIF shown');

return;
