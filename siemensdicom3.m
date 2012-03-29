function siemensdicom2(inputPath)

%clear all
%inputPath='/v/raid6/ed/MRI/3T/X100204'

path = strcat(inputPath,'/DicomData/');   % -- Location of the .IMA files -- %
path1=strcat(inputPath,'/ImageData/');
if ~exist(path1,'dir')
  [SUCCESS,MESSAGE,MESSAGEID]= mkdir(inputPath,'ImageData')  % it claims it fails, but seems to work!
end


filenotes = fopen(strcat(path1,'Notes'),'w')

names = dir(path);   % this gets list of all files in path directory

numSeries=0;
numFrames=zeros(1,10000);
numImages=zeros(1,10000); 
SeriesNumber_old=0;

for m=3:length(names)    % first two are . and ..
    m
    try f = dicominfo(strcat(path,names(m).name));
    catch disp('cant read this darn thing - keep going, may be ok')
        continue;  % want to stay in while loop, just go to next...
    end
    
    try f2 = dicomread(strcat(path,names(m).name));
    catch
        continue;
    end
    
    if numSeries < f.SeriesNumber
        numSeries = f.SeriesNumber;
    else numSeries=numSeries;
    end
    
    %%% This was for use with Siemens MRI data since the field 'AcquisitionNumber' exists
    %     try
    %     if numFrames(f.SeriesNumber) < f.AcquisitionNumber
    %         numFrames(f.SeriesNumber) = f.AcquisitionNumber;
    %     else numFrames(f.SeriesNumber)=numFrames(f.SeriesNumber);
    %     end
    %     catch
    %         continue
    %     end
    
    %%% Below is for GE data, in which the field 'AcquisitionNumber' doesn't exist
    tmp1=f.ImagesInAcquisition/f.NumberOfTemporalPositions; % added NP 041609
    %numFrames(f.SeriesNumber); % added NP 041609
    tmp2=ceil(f.InstanceNumber/tmp1); % added NP 041609
    %tmp2=f.CardiacNumberOfImages; %% NP 042909
    TMP2=1; % hmmm???, will set to 1 for now
    try
        if numFrames(f.SeriesNumber) < TMP2 % added NP 041609
            numFrames(f.SeriesNumber) = TMP2; % added NP 041609
        else numFrames(f.SeriesNumber)=numFrames(f.SeriesNumber);
        end
    catch continue
    end
    
    if numImages(f.SeriesNumber) < f.InstanceNumber
        numImages(f.SeriesNumber) = f.InstanceNumber;
    else numImages(f.SeriesNumber)=numImages(f.SeriesNumber);
    end
    
    try % this "try, catch" was added for GE data that didn't have a SeriesNumber associated with it
        %fileout =strcat(path,'series',num2strWithZeros(f.SeriesNumber),'.acq',num2strWithZeros(f.AcquisitionNumber),'.instance',num2strWithZeros(f.InstanceNumber));     %% changed since GE data always reports f.AcquisitionNumber as "1" andhere I want it to reflect each separate acquisition
        fileout = strcat(path,'series',num2strWithZeros(f.SeriesNumber),'.acq',num2strWithZeros(tmp2),'.instance',num2strWithZeros(f.InstanceNumber));
    catch
    end
    
    try copyfile(strcat(path,names(m).name), fileout);
    catch disp('not doing copyfile for'),names(m).name
        %           keyboard
    end
    
    
    if (f.InstanceNumber == 1 && f.SeriesNumber==1)   % conceivably could fail to print out this first part
        fprintf(filenotes,'Study Date: %s\t',f.StudyDate);
        fprintf(filenotes,'This file created with dicom2short %s\t',date);
        %           fprintf(filenotes,'Patient MRN: %s\t',f.PatientID);
        %           fprintf(filenotes,'Patient Name: %s\t',f.PatientsName.FamilyName);
        try fprintf(filenotes,'Pat. weight %5.1f\t\n\n',f.PatientsWeight);
        catch  disp('no weight entered')
        end
        % fprintf(filenotes,'%s\n\n',f.PatientsName.GivenName);
    end
    
end;

numFrames=numFrames(1:numSeries);
%%numSlices=numImages(1:numSeries)./numFrames; % NP for GE data these fields result in numSlices=0  % added NP 041609
tmp3=zeros(size(numImages(1:numSeries)));  % added NP 041609
tmp3(:)=f.ImagesInAcquisition/f.NumberOfTemporalPositions  % added NP 041609
numSlices=tmp3;  % added NP 041609

save numVectors numSeries numImages numFrames numSlices
% just saving if want to run part below separately...

try fclose(filenotes);
    catch('I dont know why it said invalid fidd') 
end


%---- Now start second pass through --- %

clear names
names=dir(path);
aa='';

for m=1:length(names)
   if strncmp(names(m).name,'series',6)
      aa=strvcat(aa,names(m).name);
   end
end


filelist=cellstr(aa);
filelistsorted=sort(filelist)

%dd=strmatch('sseries0007',filelist);
%ee=filelist(dd);
%clear filelist
%filelist=ee;

transposeflag=1;

for m=1:length(filelist)
       fname = char(strcat(path,filelist(m)));           
       try f = dicominfo(fname);           
          catch disp('cant read this darn thing')
                continue;  % want to stay in while loop, just go to next...
       end
       try fdata = dicomread(fname);           
          catch  continue;
       end

     if (transposeflag == 1)
           f2 = fdata';
           nRows = f.Columns;
           nColumns = f.Rows;
       else
           f2 = fdata;
       end;

% but for dynamic need 10 frames, 5 slices, then will go like 1,3,5,2,4, 1,3,...


     numPhases=1;
     %try numPhases=f.CardiacNumberOfImages % NP 041609 changed since GE scanners write out current time frame number for this field, instead of 1 frame per image--should be '1' 
     try numPhases=f.AcquisitionNumber % so it always shows 1 image per .dcm image  
        catch disp('not segemented')
      end
     if numPhases==1
        slicenum=mod(f.InstanceNumber,numSlices(f.SeriesNumber))+1;
        else
        try
          numSlices(f.SeriesNumber)=numImages(f.SeriesNumber)/numFrames(f.SeriesNumber);
          slicenum=mod(f.InstanceNumber,numSlices(f.SeriesNumber))+1;
          %%%numFrames(f.SeriesNumber)=f.CardiacNumberOfImages; % again I'll use f.AcquisitionNumber since it's always 1 with GE data... NP 041609
          numFrames(f.SeriesNumber)=f.AcquisitionNumber;
          catch slicenum=1;
        end
%     else
%       slicenum=floor((f.InstanceNumber-1)/numFrames(f.SeriesNumber))+1;
     end
    
     keyboard
     
     
     if numFrames(f.SeriesNumber)==1 && numPhases==1
       fid2 = fopen(strcat(path1,'series',num2str(f.SeriesNumber),'.slice',num2str(slicenum),'.',num2str(nRows),'x',num2str(nColumns),'x',num2str(numSlices(f.SeriesNumber))),'a','ieee-be');
       disp('writing all together if numFrames =1 for this series')
       fwrite(fid2,f2,'short');
       fclose(fid2);
     else
       if numPhases~=1
          numFrames(f.SeriesNumber)=numPhases;
       end
       fid2 = fopen(strcat(path1,'series',num2str(f.SeriesNumber),'.',num2str(nRows),'x',num2str(nColumns),'x',num2str(numFrames(f.SeriesNumber))),'a','ieee-be');
       fwrite(fid2,f2,'short');
       fclose(fid2);
     end
%keyboard

end;







function outnum=num2strWithZeros(num);
       if (num<=9)
          outnum=strcat('000',num2str(num));
       elseif ((9<num) && num<=99)
          outnum=strcat('00',num2str(num));
       elseif ((99<num) && num<=999)
          outnum=strcat('0',num2str(num));
       elseif ((999<num) && num<=9999)
          outnum=num2str(num);
       else
          disp('Problem with out of range number!!!')
       end
return