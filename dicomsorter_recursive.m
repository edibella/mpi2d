function map = dicomsorter_recursive(top_of_path2read,outputpath,seriesopt,sliceopt,map);
% DICOMSORTER is a function that sorts dicom data from a given input path
%   path2read based on their sequence/series of acquisition and slice
%   ordering. By default it only sorts on series number. Usage as follows:
%
%   dicomsorter(path2read), dicomsorter(path2read,1) - sorts on series number
%   dicomsorter(path2read, 1, 1) - sorts on series number and the
%   respective series into their slices
%
% See also SORTONSERIESNUMBER, SORTONSLICENUMBER
%
% Sathya Vijayakumar
% July 2009, UCAIR, University of Utah
% Made recursive and mkdir, EVRD 7/09

if nargin < 4
    sliceopt = 0;
end;
if nargin < 3 || (seriesopt == 0)
    seriesopt = 1;
end;
if nargin < 2 
    outputpath = pwd;
else
    if exist(outputpath)~=7
        mkdir(outputpath)
    end
end;
if nargin < 1
    errordlg('Sorry, cannot proceed without input directory');
    return;
end;
if(~exist('map'))
    map = containers.Map({'dummykey'},{zeros(1,4)}); remove(map,'dummykey');
end
    
h = dir(top_of_path2read);
debugging = 0;
h2 = waitbar(0,strcat('Please wait, reading dicom headers for directory ', top_of_path2read) );
description = '';
filenames = {};
d = {};
for ii=3:length(h)
   if (h(ii).isdir)
       map = dicomsorter_recursive( strcat(top_of_path2read,'/',h(ii).name), outputpath,map)
   else

%for i = 3:length(h)
    try
        d{end+1} = dicominfo(strcat(top_of_path2read,'/',h(ii).name));
        filenames = vertcat(filenames,strcat(top_of_path2read,'/',h(ii).name));
        
        description = [char(d{end}.SeriesDescription) '(' num2str(d{end}.SeriesNumber) ')'];
        if(~isfield(d{end},'SliceLocation'))
            d{end}.SliceLocation = 1;
        end
        if(~isKey(map,description))
            map(description) = [d{end}.SliceLocation];
        end
        list = map(description);
        index = -1;
        for i=1:length(list)
            if(list(i) == d{end}.SliceLocation)
                index = i;
                break;
            end
        end
        if(index < 1)
            newList = horzcat(map(description) ,d{end}.SliceLocation);
            newList = sort(newList);
            map(description) = newList;
        end
        
    catch ME
        if(debugging)
            disp(strcat(top_of_path2read,'/',h(end).name,'  is not a dicom file, skip and continue '))  % seems a poor way to handle this
            disp(['Would have placed it in ' description]);
        end
        % this likely only works if it is the first file listed in the
        % directory!
        %d{end}.SeriesNumber = -99;
    end
    
        
    waitbar(ii/length(h),h2);    
   end
end;
close(h2);
if(~exist('d'))
    disp(['Encounter error reading ' top_of_path2read]);
    return;
end
%if (d{1}.SeriesNumber~=-99)
   sortonseriesnumber(d,filenames,top_of_path2read,outputpath,map);
%end


