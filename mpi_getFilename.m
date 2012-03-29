function filename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions, clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF, flagTimeStamps)

if ~exist('useIntegralLinearFit')
   useIntegralLinearFit=0;
end
if ~exist('flagTimeStamps')   % a bit non-intuitive, but assume being used if not passed
   flagTimeStamps=1;
end

   filenameFirstPart=strcat(outpath,'flowvalues.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions), '.',int2str(numRadialRegions));

if fixedDelay==99 & fixedVp==99
   filenameIntermediate=filenameFirstPart;
elseif fixedDelay==99 & fixedVp~=99
   filenameIntermediate=strcat(filenameFirstPart,'_fixedVp',num2str(fixedVp));
elseif fixedDelay~=99 & fixedVp==99
   filenameIntermediate=strcat(filenameFirstPart,'_fixedDelay',num2str(fixedDelay));
elseif fixedDelay~=99 & fixedVp~=99
   filenameIntermediate=strcat(filenameFirstPart,'_fixedDelay',num2str(fixedDelay),'_fixedVp',num2str(fixedVp));
else
  disp('error in fit_useAlt writing  out parameters')
end

if exist('seriesNumAIF')
   filenameLastStage=strcat(filenameIntermediate,'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_', num2str(scaleAIF));
else
   filenameLastStage=filenameIntermediate;
end

filename=filenameLastStage;

if useIntegralLinearFit~=0    % add .linearFits.txt to all
   filename=strcat(filenameLastStage, '.linearFits');
end

if flagPixelwise==1   
   filename=strcat(filename, '.pixelwise');
end

if clusterFlag==1   
   filename=strcat(filename, '.clusters');
end

if flagTimeStamps==0   
   filename=strcat(filename, '.no_ts')
end

filename=strcat(filename, '.txt');

return
