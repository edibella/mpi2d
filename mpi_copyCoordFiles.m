function mpi_copyCoordFiles(sliceNum, studyNum)

% new take on this 11/9/04, just copy 4 files to new study and slice number

disp('WARNING!!!  THIS WILL OVERWRITE THE 4 COORD FILES, THE CONTOURS YOU HAVE DRAWN, ON THE CURRENT STUDY IF THEY EXIST.')
disp('Hit CNTRL-C to abort if you do not want them overwritten!') 

existingStudyNum=input('Enter study or series number to copy from \n')
existingSliceNum=input('Enter slice number to copy from \n')

% not logging this, figure can check which matches which if necessary by doing diffs and looking at dates...

  fileout=strcat('Output/endo_polyCoords.study',int2str(studyNum),'.slice',int2str(sliceNum));
  try copyfile(strcat('Output/endo_polyCoords.study',int2str(existingStudyNum),'.slice',int2str(existingSliceNum)), fileout);
  catch disp('cant seem to find endo Coords')
  end

  fileout=strcat('Output/epi_polyCoords.study',int2str(studyNum),'.slice',int2str(sliceNum));
  try copyfile(strcat('Output/epi_polyCoords.study',int2str(existingStudyNum),'.slice',int2str(existingSliceNum)), fileout);
  catch disp('cant seem to find epi Coords')
  end


  fileout=strcat('Output/Roi_start_angle.study',int2str(studyNum),'.slice',int2str(sliceNum));
  try copyfile(strcat('Output/Roi_start_angle.study',int2str(existingStudyNum),'.slice',int2str(existingSliceNum)), fileout);
  catch disp('cant seem to find epi Coords')
  end


  fileout=strcat('Output/blood_polyCoords.study',int2str(studyNum),'.slice',int2str(sliceNum));
  try copyfile(strcat('Output/blood_polyCoords.study',int2str(existingStudyNum),'.slice',int2str(existingSliceNum)), fileout);
  catch disp('cant seem to find epi Coords')
  end

return;
