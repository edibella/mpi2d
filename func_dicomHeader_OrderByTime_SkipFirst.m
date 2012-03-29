function func_dicomHeader_OrderByTime_SkipFirst(templatedirtxt, indata, outname, seriesnum, seriesdesc); 
%  func_dicomHeader(templatedirtxt, indata, outname, seriesnum, seriesdesc)
%   finds all files with suffix .dcm in given templatedir
% revised with more error checking  2/24/10

tmpp=strcat(templatedirtxt, '/*.dcm')  
h = dir(tmpp);
slice_loc = [];
n=0;

for i = 1:1  %42 %-- Here the assumption is for our scans, one slice per folder -- %
    h(i).name
    
    try
        d = dicominfo(strcat(templatedirtxt,'/',h(i).name));
     catch
         disp('Problem with reading dicominfo from '), d
          keyboard
     end
    
    n = n+1;
    if n == 1 || isempty(find(slice_loc == d.SliceLocation)) 
        slice_loc = [slice_loc,d.SliceLocation];
    else 
        break;
    end;
end;

slice_loc2use = slice_loc(1);   % this is like 20.582

slice_match_indices = [];
disp('Please wait... Sorting through the template file directory')
%h1 = waitbar(0,'Please wait... Sorting through the template file directory');
clear sortingValue
for i = 1:length(h)
%    waitbar(i/length(h),h1);
    d = dicominfo(strcat(templatedirtxt,'/',h(i).name));
    if d.SliceLocation == slice_loc2use
        slice_match_indices = [slice_match_indices,i];
        sortingValue(i) = str2num(d.AcquisitionTime);
	    %sortingValue(i) = d.InstanceNumber;
    else
        disp('Note that a slice location does not match!! Likely an incorrect dicom file in this folder ')
        disp('Need to resolve this')
        keyboard;
    end;
end;
%close(h1);
[sortingValue IX] = sort(sortingValue);
h = h(IX);
h = h(2:end);
slice_match_indices = slice_match_indices(IX);
slice_match_indices = slice_match_indices(2:end);

% -- Start the dicom write process -- %
if length(slice_match_indices) ~= size(indata,3)
    %errordlg('Data size (3rd dimension) mismatch, please check input dataset! Will continue','Dataset size difference');
    disp('Data size (3rd dimension) mismatch, please check input dataset! Will continue ');
end;
%h2 = waitbar(0,'Please wait... Writing out the new dicom files');
disp('Please wait... Writing out the new dicom files')

%%%% to find max intensity in the sequence
%max_int=0;
%for k = 1:size(indata,3)
%    templateFilename = strcat(templatedirtxt,'/',h(k).name);
%    tmpg=dicominfo(templateFilename);
%    
%if(double(tmpg.LargestImagePixelValue)>max_int)
%        max_int=double(tmpg.LargestImagePixelValue);
%    end    
%end;
%indata=indata*(max_int/max(indata(:)));


for k = 1:length(h)
    templateFilename = strcat(templatedirtxt,'/',h(k).name);
    
    addDicomHeader(templateFilename, indata(:,:,k), seriesnum, ...
        seriesdesc, k, outname);
   % waitbar(k/(size(indata,3)), h2);
end;
%close(h2);
return;
