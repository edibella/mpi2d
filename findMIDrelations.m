function [MIDinfolookup,SeriesNumberinfolookup] = findMIDrelations()
myprocessingDir = pwd;

try
    MIDinfolookup = containers.Map('Keytype','int32','ValueType','any');
    SeriesNumberinfolookup = containers.Map('Keytype','int32','ValueType','any');
catch
    clear temp
    temp.bla = 0;
    MIDinfolookup = containers.Map({int32(1)},{temp}); remove(MIDinfolookup,int32(1));
    SeriesNumberinfolookup = containers.Map({int32(1)},{temp}); remove(SeriesNumberinfolookup,int32(1));
end
clc
disp('-----------------------');
if(~strfind(myprocessingDir,'Processing'))
    disp('I must have the dicom folder present, as well as the RawData folder')
    return;
end
cd('..');
if(~exist('RawData'))
    cd(myprocessingDir);
    disp('I cannot continue without the RawData folder present.  Please link to it');
    return
end
if(~exist('DicomData'))
    cd(myprocessingDir);
    disp('I cannot continue without the DicomData folder present.  Please link to it or create it');
    return;
else
    cd('DicomData');
    dicomDataFolder = pwd;
    cd([myprocessingDir '/..']);
end
cd('RawData');
rawDataFolder = pwd;
rawFiles = dir('*MID*.dat');

%filter out all dat files that don't have a mat files that was created, or
%in other words, kill all dat files that haven't been converted to mat
%files
CollectedFiles = [];
OnesToKill = [];
for i=1:length(rawFiles)
    A = sscanf(rawFiles(i).name,'meas_MID%d_%s');
    if(isempty(A)), continue, end
    MID = A(1);
    sisterFiles = dir(['*MID' num2str(MID) '*']);
    hasSister = 0;
    for j=1:length(sisterFiles)
        if(~strcmp(sisterFiles(j).name,rawFiles(i).name))
            hasSister = 1;
        end
    end
    if(~hasSister)
        OnesToKill = [OnesToKill i];
    end
end

if(length(OnesToKill) == length(rawFiles))
    disp('I don''t think you''ve run any conversion script to convert the dat files to mat files');
    disp('I''ll assume that all the first MID number with it''s corresponding SeriesDescription');
    disp('gets to claim the dicoms that have that SeriesDescription');
    disp('And that all subsequent dat files will be thrown away');
else
    disp('There were dat files that didn''t get converted.  They are:');
    for i=1:length(OnesToKill)
        disp(['- ' rawFiles(OnesToKill(i)).name]);
    end
    rawFiles(OnesToKill) = [];
end




disp('-----------------------');
disp('MID: SeriesDescription');
disp('-Series Numbers found');
disp(' ');
SeriesNumbersTaken = zeros(10,1);
for i=length(rawFiles):-1:1
    cd(myprocessingDir);
    A = sscanf(rawFiles(i).name,'meas_MID%d_%s');
    if(isempty(A))
        disp(['Could not parse(' rawFiles(i).name ')']);
        disp('Please construct a struct with the following parameters: ');
        disp('mystruct.SeriesNumber = 1;');
        disp('mystruct.SliceNumber = 1;');
        disp('mystruct.SeriesDescription = ''bla'';');
        A = sscanf(rawFiles(i).name,'meas_MID%d');
        if(~isempty(A))
            disp(['MIDinfolookup(' num2str(A(1)) ') =  mystruct;']);
        else
            disp('MIDinfolookup(someMID) = mystruct;' );
        end
        continue;
    end
    clear mystruct
    MID = A(1);
    
    %extract the seriesDescripotion
    mystruct.SeriesDescription = ['' A(2:end)'];
    FIDindex = strfind(mystruct.SeriesDescription,'_FID');
    mystruct.SeriesDescription = mystruct.SeriesDescription(1:(FIDindex-1));
    description = mystruct.SeriesDescription;
    disp([num2str(MID) ' : ' mystruct.SeriesDescription]);
    
    
    cd(dicomDataFolder);
    searchString = [mystruct.SeriesDescription '(*'];
    searchString = strrep(searchString,'_','*');
    candidateDicomFolders = dir(searchString);
    if(length(candidateDicomFolders) == 0)
        disp(['      Cannot find any dicoms with Description: ' mystruct.SeriesDescription ]);
        continue;
    end
    candidateSeriesNumbers = [];
    for j=1:length(candidateDicomFolders)
        cd(candidateDicomFolders(j).name);
        dicoms = dir('*.dcm');
        dicomHeader = dicominfo(dicoms(1).name);
        cd(dicomDataFolder);
        candidateDicomFolders(j);
        dicomHeader.AcquisitionTime;
        candidateDicomFolders(j).timestamp = str2num(dicomHeader.AcquisitionTime);
        candidateDicomFolders(j).seriesNumber = dicomHeader.SeriesNumber;
        candidateSeriesNumbers = horzcat(candidateSeriesNumbers,dicomHeader.SeriesNumber);
    end
    [candidateSeriesNumbers,IX] = sort(candidateSeriesNumbers);
    candidateDicomFolders = candidateDicomFolders(IX);
    j = length(candidateDicomFolders);
    collectedSeriesNumbers = [];
    %while
    %(
    %either the next series number that is to be claimed leads to an
    %invalid index in the array that holds which series numbers have been
    %claimed 
    %   or
    %that series number has not been claimed
    %)
    %   and
    %we havent' walked off the edge of the list of series numbers to be
    %claimed (this should be the same as the first test)
    while(0 < j)      
         while(length(candidateSeriesNumbers) >= j && ...
            length(SeriesNumbersTaken) >= candidateSeriesNumbers(j) && ...
            SeriesNumbersTaken(candidateSeriesNumbers(j)) && ...
            j > 0)
            j = j - 1;
            if(j == 0)
               break;
            end
        end
        if(j == 0)
            break;
        end
        collectedSeriesNumbers = [collectedSeriesNumbers candidateSeriesNumbers(j)];
        SeriesNumbersTaken(candidateSeriesNumbers(j)) = 1;
        j = j - 1;
        %if either this is the last one or the the difference to the next
        %timestamp is more than 5 seconds we break due to that being
        %another injection
        if(j == 0 || abs(candidateDicomFolders(j).timestamp - candidateDicomFolders(j+1).timestamp) > 127)  % EVRD 1/7/11, most or all older (before June 2010
            % dicom names mean nothing (now stg 1.1 gives dicomname as instance_time). So this was 5, and could have first
            % frames in different slices different times. 
            break;
        end
    end
    
    
    mystruct.SeriesNumbers = sort(collectedSeriesNumbers);
    collectedSeriesNumbers=sort(collectedSeriesNumbers);  % EVRD 12/31/10, was flipping 1 2 3 to 3 2 1 etc. 
    disp(['      - [' num2str(mystruct.SeriesNumbers) ']']);
    MIDinfolookup(MID) = mystruct;
    clear mystruct
    mystruct.MID = MID;
    mystruct.SeriesDescription = description;
    for slice = 1:length(collectedSeriesNumbers);
        SeriesNumber = collectedSeriesNumbers(slice);
        mystruct.slice = slice;
        if(isKey(SeriesNumberinfolookup,SeriesNumber))
            disp(['Colliding data for series number: ' num2str(SeriesNumber)]);
            disp('Skipping');
            continue;
        end
        SeriesNumberinfolookup(SeriesNumber) = mystruct;
    end
end
cd(myprocessingDir);
save('MIDinfolookup.mat','MIDinfolookup','SeriesNumberinfolookup');
end