% read in raw, recon with FBP and 3 subsets, and write out as dicom

% Note, views only shared within each 72 ray time frame (could extend to do 3 other time frames for 288 rays total). 
% Script needs dicom data unpacked by dicomsorter_recursive.m    
% 
% Takes in raw (.mat) k-space file and writes out recon with appropriate dicom headers, in 
% subdirectory to the Dicom data called Recon (can be altered. Note: don't write out to same Dicom directory,
% as subsequent uses of this program could get confused)

clear

scanner='Verio'
study='P050510'
startFrame=1
addpath /Users/ed/Documents/ShareWithWindows/matlab/generalCodes



% .02 rest study
% seriesnumStart=26
% seriesnumEnd=28
% % must also change filenames
% 
% load -mat /v/raid1/jfluck/MRIdata/Cardiac/Verio/P050510/RawData/meas_MID30_CV_Radial7Off_flex_72rays3.6_FID76462_Kspace;
% MID=30
% size(kSpace)
% %    256   72 15 160 2   for P040110  and series 18
% 
% for seriesnum=seriesnumStart:seriesnumEnd
%   slice=seriesnum-seriesnumStart+1
%   kSpaceOneSlice = single(kSpace(:,:,1:3:end,startFrame:end,slice) * 1e9);  % only 3 coils! Must have been CP mode and only 3!!   
%   template_dicom_directory=strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',study,'/DicomData/CV_Radial7Off_flex_72rays3.6(',int2str(seriesnum),')')
%   output_dicom_directory=strcat('/v/raid1/ed/MRIdata/Cardiac/',scanner,'/',study,'/ReconData/CV_Radial7Off_flex_72rays3.6(' );
%   outname=strcat('MID_',int2str(MID),'_imgSW_slice_',int2str(slice),'_subset');  % will add on subset # and window size when use as filename and 
%   % used as series description
%   func_directFBPandShare_3subsets(kSpaceOneSlice, template_dicom_directory, output_dicom_directory, slice, seriesnum*1000, outname );
% % only one slice at a time for now!!!
% 
% % would like to at some point put in code to automatically run mpi2d using
% % series21 parameters...
% 
% end
   
   
   

% adeno study
seriesnumStart=29
seriesnumEnd=31
% must also change filenames

load -mat /v/raid1/jfluck/MRIdata/Cardiac/Verio/P050510/RawData/meas_MID31_CV_Radial7Off_flex_72rays5.4_ADENO_FID76463_Kspace;
MID=31

for seriesnum=seriesnumStart:seriesnumEnd
  slice=seriesnum-seriesnumStart+1
  kSpaceOneSlice = single(kSpace(:,:,1:3:end,startFrame:end,slice) * 1e9);  
  template_dicom_directory=strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',study,'/DicomData/CV_Radial7Off_flex_72rays5.4_ADENO(',int2str(seriesnum),')')
  output_dicom_directory=strcat('/v/raid1/ed/MRIdata/Cardiac/',scanner,'/',study,'/ReconData/CV_Radial7Off_flex_72rays5.4_ADENO(');
  outname=strcat('MID_',int2str(MID),'_imgSW_slice_',int2str(slice),'_subset');  % will add on subset # and window size when use as filename and 
  % used as series description
  func_directFBPandShare_3subsets(kSpaceOneSlice, template_dicom_directory, output_dicom_directory, slice, seriesnum*1000, outname );
end
   
   


% lexi study
seriesnumStart=58
seriesnumEnd=60
% must also change filenames

load -mat /v/raid1/jfluck/MRIdata/Cardiac/Verio/P050510/RawData/meas_MID75_CV_Radial7Off_flex_72rays5.4_LEXI_FID76507_Kspace;
MID=75

for seriesnum=seriesnumStart:seriesnumEnd
  slice=seriesnum-seriesnumStart+1
  kSpaceOneSlice = single(kSpace(:,:,1:3:end,startFrame:end,slice) * 1e9);  
  template_dicom_directory=strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',study,'/DicomData/CV_Radial7Off_flex_72rays5.4_LEXI(',int2str(seriesnum),')')
  output_dicom_directory=strcat('/v/raid1/ed/MRIdata/Cardiac/',scanner,'/',study,'/ReconData/CV_Radial7Off_flex_72rays5.4_LEXI(');
  outname=strcat('MID_',int2str(MID),'_imgSW_slice_',int2str(slice),'_subset');  % will add on subset # and window size when use as filename and 
  % used as series description
  func_directFBPandShare_3subsets(kSpaceOneSlice, template_dicom_directory, output_dicom_directory, slice, seriesnum*1000, outname );
end







% folder = '/v/raid1/bmatthew/MRIdata/Cardiac/Trio/P030210/Processing'
% cd(folder);
% delete('mpi2d.par');
% files = dir([folder '/*.par']);
% % for i=1:(length(files)-2)
% %     copyfile(char(files(i).name),'mpi2d.par');
% %     disp(files(i).name);
% %     mpi2d stg=[3.11]
% %     keyinput = input('press 1');
% %     while(keyinput ~= 1)
% %         pause(1);
% %     end
% % end
% 
% for i=1:length(files)
%     copyfile(char(files(i).name),'mpi2d.par');
%     mpi2d stg=[3.2 3.3 3.4 4.1 4.2]
% end
% 
% 




