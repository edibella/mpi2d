function sortonseriesnumber(d, filenames,path2read,outputpath,map)
% SORTONSERIES is a function to sort dicom images based on series and write
%   out files into separate folders. The inputs are the dicom headers d,
%   the input directory path2read and the directory where the sorted files
%   are to be copied outputpath.
%   
% See also SORTONSLICENUMBER, DICOMSORTER

% Sathya Vijayakumar
% July 2009 UCAIR, University of Utah

if nargin < 3
    outputpath = pwd;
end; 
if nargin < 2
    errordlg('Sorry, cannot proceed without input data directory');
    return
end;
if nargin < 1
    errordlg ('Sorry, cannot proceed without dicome headers');
    return;
end;

h = filenames;
seriesnum = [];
k = 0;


j2 = waitbar(0,'Please wait, copying files to respective folders');
for i = 1:length(d)
    if(~isfield(d{i},'SeriesNumber')  || d{i}.SeriesNumber==-99)
        continue;
    end
    
    %find the slice number
    %if there are more than 1 slice in this injection, multiply the series number by 10 and add the slice number
    description = [char(d{end}.SeriesDescription) '(' num2str(d{end}.SeriesNumber) ')'];
    if(~isKey(map,description))
        disp('Problem here.  Couldn''t find the description after I already should have constructed the slice location map');
        keyboard
    end
    list = map(description);
    changedSomething = 0;
    if(length(list) > 1)
        sliceNumber = find(list == d{i}.SliceLocation);
        if(isempty(sliceNumber))
            keyboard
        end
        d{i}.SeriesNumber = 10*d{i}.SeriesNumber + sliceNumber;
        changedSomething = 1;
    end
    
    folderDescription = [d{i}.SeriesDescription '(' num2str(d{i}.SeriesNumber) ')'];
    if(exist([outputpath '/' folderDescription]) ~= 7)
        mkdir([outputpath '/' folderDescription]);
    end
    
    %new file name
    [instance, errmsg] = sprintf('%04d',d{i}.InstanceNumber);
    if(~isfield(d{i},'AcquisitionTime'))
        time = num2str(i);
    else
        time = d{i}.AcquisitionTime;
    end
    newFileName = [instance ' - ' time];
    if(isempty(strfind(newFileName,'.dcm')))
        newFileName = [newFileName '.dcm'];
    end
    
    if(changedSomething)
        dicomwrite(dicomread(filenames{i}),[outputpath '/' folderDescription '/' newFileName],d{i});
    else
        copyfile(filenames{i},[outputpath '/' folderDescription '/' newFileName]);
    end
    
    
    
    
%     
%     
%     if isempty(find(seriesnum == d{i}.SeriesNumber))
%         k = k+1;
%         seriesnum(k) = d{i}.SeriesNumber;
%         seriesname{k} = d{i}.SeriesDescription;
%         
%         description = [char(d{i}.SeriesDescription) '(' num2str(d{i}.SeriesNumber) ')'];
%         if(~isKey(map,description))
%             map(description) = {d{i}.SliceLocation};
%         end
%         list = map(description);
%         index = -1;
%         for ii=1:length(list)
%             if(list(ii) == d{i}.SliceLocation)
%                 index = ii;
%                 break;
%             end
%         end
%         if(index < 1)
%             map(description) = horzcat(map(description) ,d{i}.SliceLocation);
%             index = length(map(description));
%         end
%         if(length(map(description)) > 1)
%             dirname{k} = strcat(outputpath,'/',seriesname{k},'(',num2str(seriesnum(k)+index-1),')');
%         else
%             dirname{k} = strcat(outputpath,'/',seriesname{k},'(',num2str(seriesnum(k)),')');
%         end
%         if(~isdir(dirname{k}))
%             mkdir(dirname{k});
%         end
% %         [SUCCESS,message,MESSAGEID] = mkdir(dirname{k})       %troubleshooting
% %         if ~isempty(message)
% %             dirname{k}
% %             keyboard
% %         end
%         k1 = k;
%     else 
%         k1 = find(seriesnum == d{i}.SeriesNumber);
%     end;
%     num = 1;
%     try
%         targetFileName = d{i}.AcquisitionTime;
%     catch ME
%         targetFileName = h(i+2).name;
%     end
%     if(exist(strcat(dirname{k1},'/',targetFileName,'.dcm')))
%         while(exist(strcat(dirname{k1},'/',targetFileName,num2str(num),'.dcm')))
%             num = num+1;
%         end
%         targetFileName = strcat(dirname{k1},'/',targetFileName,num2str(num),'.dcm');
%     else
%         targetFileName = strcat(dirname{k1},'/',targetFileName,'.dcm');
%     end
%     copyfile(strcat(path2read,'/',h(i+2).name), targetFileName);
%     
    waitbar(i/length(d),j2);
end

close(j2);

        