function siemensdicom1(inputPath, ouputPath)

path = strcat(inputPath,'/DicomData/');   % -- Location of the .dcm files -- %
path1=strcat(ouputPath,'/ImageData/');   % -- Location of the .IMA files -- %
names = dir([inputPath '/*.dcm']);
if length(names) > 0
    path = [inputPath '/'];
end

if ~exist(path1,'dir')
  [SUCCESS,MESSAGE,MESSAGEID]= mkdir(inputPath,'ImageData');  % it claims it fails, but seems to work!
end


filenotes = fopen(strcat(path1,'Notes'),'w');

names = dir([path '/*.dcm']);   % this gets list of all files in path directory

numSeries=0; numFrames=zeros(1,10000);
numImages=zeros(1,10000); 

for m=1:length(names)    % first two are . and ..
    m
    try f = dicominfo(strcat(path,names(m).name));
    catch
        disp('cant read this darn thing - keep going, may be ok')
        continue;  % want to stay in while loop, just go to next...
    end
%     try f2 = dicomread(strcat(path,names(m).name));           
%     catch  continue;
%     end
    
    numSeries = max(f.SeriesNumber,numSeries);

    try
       numFrames(f.SeriesNumber) = max(f.AcquisitionNumber,numFrames(f.SeriesNumber));
    catch continue
    end

    numImages(f.SeriesNumber) = max(numImages(f.SeriesNumber),f.InstanceNumber);
   

end;

numFrames=numFrames(1:numSeries);
numSlices=numImages(1:numSeries)./numFrames;

save numVectors numSeries numImages numFrames numSlices
% just saving if want to run part below separately...

try fclose(filenotes);
    catch('I dont know why it said invalid fidd') 
end







%---- Now start second pass through --- %

clear names
filelist=dir([path '/*.dcm']);

transposeflag=1;

for m=1:length(filelist)
    fname = [path filelist(m).name];          
    try f = dicominfo(fname);           
    catch
        disp(['Cannot read file:' filelist(m).name])
        continue;  % want to stay in while loop, just go to next...
    end


    try fdata = dicomread(fname); catch  continue;  end

    if (transposeflag == 1)
        f2 = fdata';
        nRows = f.Columns;
        nColumns = f.Rows;
    else
        f2 = fdata;
    end;

% but for dynamic need 10 frames, 5 slices, then will go like 1,3,5,2,4, 1,3,...


     numPhases=1;
     try numPhases=f.CardiacNumberOfImages; catch disp('not segemented'); end
     if numPhases==1
         slicenum=mod(f.InstanceNumber,numSlices(f.SeriesNumber))+1;
     else
        try
            numSlices(f.SeriesNumber)=numImages(f.SeriesNumber)/numFrames(f.SeriesNumber);
            slicenum=mod(f.InstanceNumber,numSlices(f.SeriesNumber))+1;
            numFrames(f.SeriesNumber)=f.CardiacNumberOfImages;
        catch
            slicenum=1;
        end
     end

     if numFrames(f.SeriesNumber)==1 && numPhases==1
        fid2 = fopen(strcat(path1,'series',num2str(f.SeriesNumber),'.',num2str(nRows),'x',num2str(nColumns),'x',num2str(numSlices(f.SeriesNumber))),'a','ieee-be');
        disp('writing all together if numFrames =1 for this series')
        fwrite(fid2,f2,'short');
        fclose(fid2);
     else
        if numPhases~=1
            numFrames(f.SeriesNumber)=numPhases;
        end
        fid2 = fopen(strcat(path1,'series',num2str(f.SeriesNumber),'.slice',num2str(slicenum),'.',num2str(nRows),'x',num2str(nColumns),'x',num2str(numFrames(f.SeriesNumber))),'a','ieee-be');
        if fid2 < 0
            stoppage = 0;
        end
        fwrite(fid2,f2,'short');
        fclose(fid2);
     end


end;