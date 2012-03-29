function sortDicoms(HeadOfDicoms, TargetFolder)

%gather all the possible dicom files
mypwd = pwd;
cd(HeadOfDicoms)
possibleFiles = gatherFilesRecusivly();
cd(mypwd);
%find all the slice locations and which files are dicoms
temp = containers.Map({0},{{''}});
sliceLocations = containers.Map({''},{temp}); remove(sliceLocations,'');
SeriesNumbersObserved = [];
waith = waitbar(0,'Finding Slice Locations and sorting');
Header = dicominfo([HeadOfDicoms '/' possibleFiles{ceil(rand(1,1)*length(possibleFiles))}]);
Headers = containers.Map({possibleFiles{1}},{Header}); 
for filei=1:length(possibleFiles)
    try
        Header = dicominfo([HeadOfDicoms '/' possibleFiles{filei}]);
    catch ME
        disp([possibleFiles{filei} ' was not a dicom file']);
        continue;
    end
    Headers(possibleFiles{filei}) = Header;
    if(~isfield(Header,'SeriesNumber'))
        Header.SeriesNumber = 1;
    end
    if(~isfield(Header,'SeriesDescription'))
        continue;
    end
    SeriesNumbersObserved(Header.SeriesNumber) = 1;
    %see if this description(folder name) has been observed
    description = [Header.SeriesDescription '(' num2str(Header.SeriesNumber) ')'];
    if(~isKey(sliceLocations,description))   %make an entry
        temp = containers.Map({0},{{''}}); remove(temp,0);
        sliceLocations(description) = temp;
    end
    %check to see if this location has been observed
    DicomFileNameByLocation = sliceLocations(description);
    if(~isfield(Header,'SliceLocation'))  %make sure there's a location.  Those without a SliceLocation field should be grouped together
        Header.SliceLocation = 1;
    end
    if(~isKey(DicomFileNameByLocation,Header.SliceLocation))  %if we haven't observed this location before add an entry
        DicomFileNameByLocation(Header.SliceLocation) = {};
    end
    %add this file by it's slice location in the map which holds those with that description
    FileNames = DicomFileNameByLocation(Header.SliceLocation);
    DicomFileNameByLocation(Header.SliceLocation) = vertcat(DicomFileNameByLocation(Header.SliceLocation),possibleFiles{filei});
    sliceLocations(description) = DicomFileNameByLocation;
    waitbar(filei/(length(possibleFiles)),waith,[strrep(description,'_','\_') '[' num2str(round(cell2mat(keys(DicomFileNameByLocation)))) ']']);
end

folders = keys(sliceLocations);
%condensing

%keyboard  %EVRD 3/2012
% looks like folders is alphabetical:
%folders(:)

% aa=(values(DicomFileNameByLocation))
% bb=cat(1,aa{:})
% cell2mat(bb)
% 
% dd=(values(sliceLocations(cell2mat(folders(1)))))
% bb=cat(1,dd{:})
% cell2mat(bb)


for folderi = 1:length(folders)
    description = folders{folderi};
    
    DicomFileNameByLocation = sliceLocations(description);
    Locations = sort(cell2mat(keys(DicomFileNameByLocation)));
    if(length(Locations) > 1)
        for slice = 2:length(Locations)
            dicomFiles = DicomFileNameByLocation(Locations(slice));
            if(length(dicomFiles) == 1)
                DicomFileNameByLocation(Locations(1)) = vertcat(DicomFileNameByLocation(Locations(1)),dicomFiles);
                remove(DicomFileNameByLocation,Locations(slice));
                sliceLocations(description) = DicomFileNameByLocation;
            end
        end
        % seems like the above takes a folder with multiple slices, and if
        % there are no time frames, combines it into one sliceLocations
        % entry (?)  EVRD 3/12
        
    end
end


filesCopied = 0;
for folderi = 1:length(folders)
    description = folders{folderi};
    waitbar(filesCopied/length(possibleFiles),waith,['Writing files to ' strrep(description,'_','\_')]);
    DicomFileNameByLocation = sliceLocations(description);
    Locations = sort(cell2mat(keys(DicomFileNameByLocation)));
    %all dicom files here are implied to have the same series number due to the description having the series number in it
    changedSomething = 0; clear NewSeriesNumber
%     if(length(Locations) > 1)  % why give multi-slcie dynamic new series number at the end??  EVRD
%         changedSomething = 1;
%         %now check to see if we can take up length(Locations) of contigous series numbers
%         %length(SeriesNumbersObserved) == highest series number observed
%         NewSeriesNumber = length(SeriesNumbersObserved)+1;
%         SeriesNumbersObserved(NewSeriesNumber:(NewSeriesNumber+length(Locations)-1)) = 1;
%     end
    %keyboard
    for slice = 1:length(Locations)
        dicomFiles = DicomFileNameByLocation(Locations(slice));
        for dicomi = 1:length(dicomFiles)
            Header = Headers(dicomFiles{dicomi});
            [instance, errmsg] = sprintf('%04d',Header.InstanceNumber);
            if(~isfield(Header,'AcquisitionTime'))
                time = num2str(slice);
            else
                time = Header.AcquisitionTime;
            end
            newFileName = [instance '_' time '.dcm'];
            if(changedSomething && exist('NewSeriesNumber'))
                newFolderName = strrep(description,['(' num2str(Header.SeriesNumber) ')'],['(' num2str(NewSeriesNumber+slice-1) ')']);
                if(~exist([TargetFolder '/' newFolderName]))
                    mkdir([TargetFolder '/' newFolderName]);
                end
                newFileName = [newFolderName '/' newFileName];
                Header.SeriesNumber = NewSeriesNumber+slice-1;
                image = dicomread([HeadOfDicoms '/' dicomFiles{dicomi}]);
				try
	                dicomwrite(image,[TargetFolder '/' newFileName],Header);
				catch ME
					disp(['Tried to write dicom(' TargetFolder '/' newFileName ') but failed']);
					continue;
				end
            else
                if(~exist([TargetFolder '/' description]))
                    mkdir([TargetFolder '/' description]);
                end
                newFileName = [description '/' newFileName];
                copyfile([HeadOfDicoms '/' dicomFiles{dicomi}],[TargetFolder '/' newFileName],'f');
            end
        end
        filesCopied = filesCopied+length(dicomFiles);
        if(changedSomething)
            disp(['Copied ' num2str(length(dicomFiles)) ' into ' newFileName(1:(end-25))]);
        else
            disp(['Copied ' num2str(length(dicomFiles)) ' into ' description]);
        end
    end
end
close(waith);

