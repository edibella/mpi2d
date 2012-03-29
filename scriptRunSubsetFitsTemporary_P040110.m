% read in raw, recon with FBP and 3 subsets, and write out as dicom

% Note, views only shared within each 72 ray time frame (could extend to do 3 other time frames for 288 rays total). 
% Script needs dicom data unpacked by dicomsorter_recursive.m    
% 
% Takes in raw (.mat) k-space file and writes out recon with appropriate dicom headers, in 
% subdirectory to the Dicom data called Recon (can be altered. Note: don't write out to same Dicom directory,
% as subsequent uses of this program could get confused)

clear

scanner='Verio'
Patient='P040110'
startFrame=1
%addpath /Users/ed/Documents/ShareWithWindows/matlab/generalCodes
addpath /v/raid1/ed/src/matlab/Radial
addpath /v/raid1/ed/src/matlab/generalCodes
addpath /v/raid1/bmatthew/Code
%addpath('/v/raid1/npack/Code/mpi2dCode_THK_092209/')

doReconFlag=0  % in case recons already done, don't need to repeat. 
% what if mother not analyzed?
doProcessing=1  % only do this if "mother" series have been registered and segmented. 

% first injection, so we can get T1_0
seriesnumStart=14
seriesnumEnd=15
% must also change filenames
scan_name='CV_Radial7Off_flex_72rays1_5(';
% 
if doReconFlag==1
    load -mat '/v/raid1/bmatthew/MRIdata/Cardiac/Verio/P040110/RawData/meas_MID103_CV_Radial7Off_flex_72rays1_5_FID68831_Kspace.mat';
    MID=103
end
% size(kSpace)
%    256   72 15 160 2   for P040110  and series 18

% for seriesnum=seriesnumStart:seriesnumEnd
%   slice=seriesnum-seriesnumStart+1
%   if doReconFlag
%   kSpaceOneSlice = single(kSpace(:,:,1:3:end,startFrame:end,slice) * 1e9);    
%   template_dicom_directory=strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',Patient,'/DicomData/CV_Radial7Off_flex_72rays1_5(',int2str(seriesnum),')')
%   output_dicom_directory=strcat('/v/raid1/ed/MRIdata/Cardiac/',scanner,'/',Patient,'/ReconData/CV_Radial7Off_flex_72rays1_5(' );
%   outname=strcat('MID_',int2str(MID),'_imgSW_slice_',int2str(slice),'_subset');  % will add on subset # and window size when use as filename and 
%   % used as series description
%   func_directFBPandShare_3subsets(kSpaceOneSlice, template_dicom_directory, output_dicom_directory, slice, seriesnum*1000, outname );
% % only one slice at a time for now!!!
%   end
% 
% if doProcessing
%     % get current directory so we can return to it afterwards
%     cd(strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',Patient,'/Processing'));
%     for subset=1:3
%         folder=strcat(scan_name,int2str(seriesnum*1000+subset),')');
%         %mpi2d stg=0.3 DoThisFolder='CV_Radial7Off_flex_72rays4.4(18002)'
%         eval(['mpi2d stg=0.3 DoThisFolder=''' folder ''''])
%         delete('mpi2d.par')
%         parfile=strcat('rad.series',int2str(seriesnum*1000+subset),'.slice',int2str(slice),'.subset',int2str(subset),'.par')
%         copyfile(parfile,'mpi2d.par')
%         %copyfile('rad.series18002.slice1.subset2.par','mpi2d.par')
%         mpi2d stg=[3.2 3.3]
%     end
%     precontrastT1=fit_subsets(Patient, seriesnum*1000, slice);
% end
% 
% end

% to pull out precontrastT1 for other studies:
T10=[];
for seriesnum=seriesnumStart:seriesnumEnd
  slice=seriesnum-seriesnumStart+1
load(strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',Patient,'/Processing/Output/MultiSRT_Vals_M0_fA_T1_chi2err.series',int2str(seriesnum*1000),'.slice',int2str(slice),'.mat'));
T10=[T10 Params.EstPrecontrast];
end
precontrastT1=median(T10)


% % .02 rest study
% seriesnumStart=18
% seriesnumEnd=19
% % must also change filenames
% scan_name='CV_Radial7Off_flex_72rays4.4('
% 
% if doReconFlag==1
%    load -mat '/v/raid1/bmatthew/MRIdata/Cardiac/Verio/P040110/RawData/meas_MID105_CV_Radial7Off_flex_72rays4.4_FID68833_Kspace';
%    MID=105
% end
% 
% for seriesnum=seriesnumStart:seriesnumEnd
%   slice=seriesnum-seriesnumStart+1
%   if doReconFlag
%   kSpaceOneSlice = single(kSpace(:,:,1:3:end,startFrame:end,slice) * 1e9);  % only 3 coils! Must have been CP mode and only 3!!   
%   template_dicom_directory=strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',Patient,'/DicomData/CV_Radial7Off_flex_72rays4.4(',int2str(seriesnum),')')
%   output_dicom_directory=strcat('/v/raid1/ed/MRIdata/Cardiac/',scanner,'/',Patient,'/ReconData/CV_Radial7Off_flex_72rays4.4(' );
%   outname=strcat('MID_',int2str(MID),'_imgSW_slice_',int2str(slice),'_subset');  % will add on subset # and window size when use as filename and 
%   % used as series description
%   func_directFBPandShare_3subsets(kSpaceOneSlice, template_dicom_directory, output_dicom_directory, slice, seriesnum*1000, outname );
% % only one slice at a time for now!!!
%    end
% 
% if doProcessing
%     % get current directory so we can return to it afterwards
%     cd(strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',Patient,'/Processing'));
%     for subset=1:3
%         folder=strcat(scan_name,int2str(seriesnum*1000+subset),')');
%         %mpi2d stg=0.3 DoThisFolder='CV_Radial7Off_flex_72rays4.4(18002)'
%         eval(['mpi2d stg=0.3 DoThisFolder=''' folder ''''])
%         delete('mpi2d.par')   %needed because copyfile won't overwrite if file has a different owner (even if e.g. group writable)
%         parfile=strcat('rad.series',int2str(seriesnum*1000+subset),'.slice',int2str(slice),'.subset',int2str(subset),'.par')
%         copyfile(parfile,'mpi2d.par')
%         %copyfile('rad.series18002.slice1.subset2.par','mpi2d.par')
%         mpi2d stg=[3.2 3.3]
%     end
%     fit_subsets(Patient, seriesnum*1000, slice, precontrastT1);
% end
% 
% end
   
   

% adeno study
seriesnumStart=22
seriesnumEnd=23
% must also change filenames
scan_name='CV_Radial7Off_flex_72rays6.5('

if doReconFlag==1
load -mat /v/raid1/gadluru/MRIdata/Cardiac/Verio/P040110/RawData/meas_MID107_CV_Radial7Off_flex_72rays6.5_FID68835_Kspace;
MID=107
end

for seriesnum=seriesnumStart:seriesnumEnd
  slice=seriesnum-seriesnumStart+1
 if doReconFlag
  kSpaceOneSlice = single(kSpace(:,:,1:3:end,startFrame:end,slice) * 1e9);  
  template_dicom_directory=strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',Patient,'/DicomData/CV_Radial7Off_flex_72rays6.5(',int2str(seriesnum),')')
  output_dicom_directory=strcat('/v/raid1/ed/MRIdata/Cardiac/',scanner,'/',Patient,'/ReconData/CV_Radial7Off_flex_72rays6.5(');
  outname=strcat('MID_',int2str(MID),'_imgSW_slice_',int2str(slice),'_subset');  % will add on subset # and window size when use as filename and 
  % used as series description
  func_directFBPandShare_3subsets(kSpaceOneSlice, template_dicom_directory, output_dicom_directory, slice, seriesnum*1000, outname );
     end

if doProcessing
    % get current directory so we can return to it afterwards
    cd(strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',Patient,'/Processing'));
    for subset=1:3
        folder=strcat(scan_name,int2str(seriesnum*1000+subset),')');
        %mpi2d stg=0.3 DoThisFolder='CV_Radial7Off_flex_72rays4.4(18002)'
        eval(['mpi2d stg=0.3 DoThisFolder=''' folder ''''])
        delete('mpi2d.par')   %needed because copyfile won't overwrite if file has a different owner (even if e.g. group writable)
        parfile=strcat('rad.series',int2str(seriesnum*1000+subset),'.slice',int2str(slice),'.subset',int2str(subset),'.par')
        copyfile(parfile,'mpi2d.par')
        %copyfile('rad.series18002.slice1.subset2.par','mpi2d.par')
        mpi2d stg=[3.2 3.3]
    end
    fit_subsets(Patient, seriesnum*1000, slice, precontrastT1);
end
end
   
   


% lexi study
seriesnumStart=65
seriesnumEnd=66
% must also change filenames
scan_name='CV_Radial7Off_flex_72rays6.5LEXI('

if doReconFlag==1
    load -mat /v/raid1/gadluru/MRIdata/Cardiac/Verio/P040110/RawData/meas_MID134_CV_Radial7Off_flex_72rays6.5LEXI_FID68862_Kspace;
    MID=134
end

for seriesnum=seriesnumStart:seriesnumEnd
  slice=seriesnum-seriesnumStart+1
  if doReconFlag==1
  kSpaceOneSlice = single(kSpace(:,:,1:3:end,startFrame:end,slice) * 1e9);  
  template_dicom_directory=strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',Patient,'/DicomData/CV_Radial7Off_flex_72rays6.5LEXI(',int2str(seriesnum),')')
  output_dicom_directory=strcat('/v/raid1/ed/MRIdata/Cardiac/',scanner,'/',Patient,'/ReconData/CV_Radial7Off_flex_72rays6.5LEXI(');
  outname=strcat('MID_',int2str(MID),'_imgSW_slice_',int2str(slice),'_subset');  % will add on subset # and window size when use as filename and 
  % used as series description
  func_directFBPandShare_3subsets(kSpaceOneSlice, template_dicom_directory, output_dicom_directory, slice, seriesnum*1000, outname );
  end
  if doProcessing
    % get current directory so we can return to it afterwards
    cd(strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',Patient,'/Processing'));
    for subset=1:3
        folder=strcat(scan_name,int2str(seriesnum*1000+subset),')');
        %mpi2d stg=0.3 DoThisFolder='CV_Radial7Off_flex_72rays4.4(18002)'
        eval(['mpi2d stg=0.3 DoThisFolder=''' folder ''''])
        delete('mpi2d.par')   %needed because copyfile won't overwrite if file has a different owner (even if e.g. group writable)
        parfile=strcat('rad.series',int2str(seriesnum*1000+subset),'.slice',int2str(slice),'.subset',int2str(subset),'.par')
        copyfile(parfile,'mpi2d.par')
        %copyfile('rad.series18002.slice1.subset2.par','mpi2d.par')
        mpi2d stg=[3.2 3.3]
    end
    fit_subsets(Patient, seriesnum*1000, slice, precontrastT1);
end
end




% lexi2 study
seriesnumStart=67
seriesnumEnd=68
% must also change filenames
scan_name='CV_Radial7Off_flex_72rays6.5LEXI2('

if doReconFlag
load -mat /v/raid1/gadluru/MRIdata/Cardiac/Verio/P040110/RawData/meas_MID135_CV_Radial7Off_flex_72rays6.5LEXI2_FID68863_Kspace;
MID=135
end

for seriesnum=seriesnumStart:seriesnumEnd
  slice=seriesnum-seriesnumStart+1
  if doReconFlag
  kSpaceOneSlice = single(kSpace(:,:,1:3:end,startFrame:end,slice) * 1e9);  
  template_dicom_directory=strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',Patient,'/DicomData/CV_Radial7Off_flex_72rays6.5LEXI2(',int2str(seriesnum),')')
  output_dicom_directory=strcat('/v/raid1/ed/MRIdata/Cardiac/',scanner,'/',Patient,'/ReconData/CV_Radial7Off_flex_72rays6.5LEXI2(');
  outname=strcat('MID_',int2str(MID),'_imgSW_slice_',int2str(slice),'_subset');  % will add on subset # and window size when use as filename and 
  % used as series description
  func_directFBPandShare_3subsets(kSpaceOneSlice, template_dicom_directory, output_dicom_directory, slice, seriesnum*1000, outname );
  end
  if doProcessing
    % get current directory so we can return to it afterwards
    cd(strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',Patient,'/Processing'));
    for subset=1:3
        folder=strcat(scan_name,int2str(seriesnum*1000+subset),')');
        %mpi2d stg=0.3 DoThisFolder='CV_Radial7Off_flex_72rays4.4(18002)'
        eval(['mpi2d stg=0.3 DoThisFolder=''' folder ''''])
        delete('mpi2d.par')   %needed because copyfile won't overwrite if file has a different owner (even if e.g. group writable)
        parfile=strcat('rad.series',int2str(seriesnum*1000+subset),'.slice',int2str(slice),'.subset',int2str(subset),'.par')
        copyfile(parfile,'mpi2d.par')
        %copyfile('rad.series18002.slice1.subset2.par','mpi2d.par')
        mpi2d stg=[3.2 3.3]
    end
    fit_subsets(Patient, seriesnum*1000, slice, precontrastT1);
  end
end





% folder = '/v/raid1/bmatthew/MRIdata/Cardiac/Trio/P040110/Processing'
% cd(folder);
% delete('mpi2d.par');
%mpi2d stg=0.2

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




