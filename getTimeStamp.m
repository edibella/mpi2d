function getTimeStamp(Patient,Scanner,Series,nSlices,nFrames)
 %%% ie. getTimeStamp('P080609','Trio',16,1,90)

%%% This routine reads in dicom file headers to determine the
%%% AcquisitionTime stamps for each image in a perfusion series.    
%%%                                                                     
%%% Before running, the user must define the following fields:



for SLICE=1:nSlices
    timeStamp=zeros(nFrames,1);
    i=1:nFrames; ii=SLICE:nSlices:nFrames*nSlices;
    acq=sprintf('%.4d',i); inst=sprintf('%.4d',ii);  % this attaches zeros to the front of the integer numbers, i, to make them have 4 digits
    counter=0;
    for n=1:4:length(acq)
        counter=counter+1;
        if Series >= 10 % ie. greater than 10 AND less than 100----(If greater than 100, remove another zero from the "series00" line below)
            D=dicominfo(strcat('/v/raid1/npack/MRIdata/Cardiac/',Scanner,'/',Patient,'/DicomData/series00',int2str(Series),'.acq',acq(n),acq(n+1),acq(n+2),acq(n+3),'.instance',inst(n),inst(n+1),inst(n+2),inst(n+3)));
        else
            D=dicominfo(strcat('/v/raid1/npack/MRIdata/Cardiac/',Scanner,'/',Patient,'/DicomData/series000',int2str(Series),'.acq',acq(n),acq(n+1),acq(n+2),acq(n+3),'.instance',inst(n),inst(n+1),inst(n+2),inst(n+3)));
        end
        tmp=D.AcquisitionTime; timeStamp(counter)=str2num(tmp);
    end
    r=find(timeStamp(:)>.9*max(timeStamp(:)));
    timeStamp=timeStamp(r);

%%%% to see the original timeStamps with erroneous shifts
% figure,plot(timeStamp(:),'b-o')
%%%% Next to correct for any erroneous time shifts...
    for i=1:size(timeStamp,1)-1
        %%% Below the threshold of total time gaps is set to 10x the mean.  However, if one jump is so big that it brings the mean to a
        %%% very high level, the other large (but comparatively smaller jumps will not be corrected.  To solve this, try changing the
        %%% threshod to 1*mean(...) below.
        if timeStamp(i+1)-timeStamp(i) > 10*mean(timeStamp(2:end)-timeStamp(1:end-1)) % the threshold value of '10' is arbitrary chosen and can be changed
           timeStamp(1:i)=timeStamp(1:i) + (timeStamp(i+1)-timeStamp(i)-0.8); % 0.8 is an average time jump between points
        end
    end
    tmp2(:,SLICE)=timeStamp;
end

timeStamp=tmp2-min(min(tmp2));
timeStampoutfile=strcat('/v/raid1/npack/Processing/',Patient,'/timeStampSer',int2str(Series),'.mat')
save(timeStampoutfile, 'timeStamp')

MeanDELTA_T=mean(timeStamp(2:end,1)-timeStamp(1:end-1,1))  % I'm just choosing the timeStamp from slice 1

%%%% to see the corrected and zeroed timeStamps compared to the original ones...
% hold on, plot(timeStamp,'g'),xlabel('Image Frames'),ylabel('Time (sec)')
% legend('Shifted Time Stamp','Corrected Time Stamp')
