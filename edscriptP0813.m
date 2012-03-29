% read in raw, recon with FBP and 3 subsets, and write out as dicom

% Note, views only shared within each 72 ray time frame (could extend to do 3 other time frames for 288 rays total). 
% Script needs dicom data unpacked by dicomsorter_recursive.m    
% 
% Takes in raw (.mat) k-space file and writes out recon with appropriate dicom headers, in 
% subdirectory to the Dicom data called Recon (can be altered. Note: don't write out to same Dicom directory,
% as subsequent uses of this program could get confused)

clear

scanner='Trio'
Patient='P081309'
startFrame=1
%addpath /Users/ed/Documents/ShareWithWindows/matlab/generalCodes
addpath /v/raid1/ed/src/matlab/Radial
addpath /v/raid1/ed/src/matlab/generalCodes
addpath /v/raid1/bmatthew/Code
RawOwner = 'bmatthew';
%addpath('/v/raid1/npack/Code/mpi2dCode_THK_092209/')
flipFlag = 1;
doReconFlag=1;  % in case recons already done, don't need to repeat. 
% what if mother not analyzed?
% do all of these - writes to my directory. 
% then cp -pr
% /v/raid1/ed/MRIdata/Cardiac/Verio/P041910/ReconData/*00[1-3]*
% /v/raid1/gadluru/MRIdata/Cardiac/Verio/P041910/ReconData/Subsets
doProcessing=1  % only do this if "mother" series have been registered and segmented. 
modifyWindowSize = 0.8;
% first injection, so we can get T1_0
seriesnumStart=18
seriesnumEnd=20
% % must also change filenames
scan_name='CV_Radial7Off_1_5(';
MID=34;

if doReconFlag==1
    load(['/v/raid1/' RawOwner '/MRIdata/Cardiac/' scanner '/' Patient '/RawData/meas_MID34_CV_Radial7Off_1_5_FID6338_Kspace.mat'],'-mat');
end
% size(kSpace)
%    256   72 15 160 2   for P040110  and series 18
% 
for seriesnum=seriesnumStart:seriesnumEnd
	slice=seriesnum-seriesnumStart+1
    template_dicom_directory=strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',Patient,'/DicomData/',scan_name,int2str(seriesnum),')')
    output_dicom_directory=strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',Patient,'/ReconData/Subsets/',scan_name );
    outname=strcat('MID_',int2str(MID),'_imgSW_slice_',int2str(slice),'_subset');  % will add on subset # and window size when use as filename and 

    if doReconFlag
        kSpaceOneSlice = single(kSpace(:,:,1:3:end,startFrame:end,slice) * 1e9);    
        % used as series description
        func_directFBPandShare_3subsets(kSpaceOneSlice, template_dicom_directory, output_dicom_directory, slice, seriesnum*1000, outname );
        % only one slice at a time for now!!!
    end
  
    if flipFlag
        for subset=1:3
            % get directory name
            dirname=strcat(output_dicom_directory,int2str(seriesnum*1000+subset),')/')
            files=dir(strcat(dirname,'*.dcm'))
            for bb=1:length(files) 
                tmpp=dicomread( strcat(dirname,files(bb).name));
                %tmpp=flipud(tmpp); 
                %imgrr(:,:,bb)=rot90(tmpp);
                imgrr(:,:,bb)=flipud(tmpp);
            end

            seriesdesc = outname;
            seriesdesc_out = strcat(seriesdesc,int2str(subset),'_window',int2str(10*modifyWindowSize));
            % write out as dicom:
            func_dicomHeader_OrderByTime_SkipFirst(template_dicom_directory, imgrr, [output_dicom_directory num2str(seriesnum*1000+subset) ')/' outname num2str(subset) '_window' int2str(10*modifyWindowSize) '_'], 1000*seriesnum+subset, seriesdesc_out);
        end
    end
    %keyboard
    if doProcessing
        % get current directory so we can return to it afterwards
        cd(strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',Patient,'/Processing'));
        for subset=1:3
            folder=strcat(scan_name,int2str(seriesnum*1000+subset),')');
            %mpi2d stg=0.3 DoThisFolder='CV_Radial7Off_flex_72rays4.4(18002)'
            eval(['mpi2d stg=0.3 DoThisFolder=''' folder ''''])
            delete('mpi2d.par')
            parfile=strcat('rad.series',int2str(seriesnum*1000+subset),'.slice',int2str(slice),'.subset',int2str(subset),'.par')
            copyfile(parfile,'mpi2d.par')
            %copyfile('rad.series18002.slice1.subset2.par','mpi2d.par')
            mpi2d stg=[3.2 3.3]
        end
        precontrastT1=fit_subsets(scanner,Patient, seriesnum*1000, slice);
    end

end

%keyboard

% to pull out precontrastT1 for other studies:
T10=[];
for seriesnum=seriesnumStart:seriesnumEnd
  slice=seriesnum-seriesnumStart+1
load(strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',Patient,'/Processing/Output/MultiSRT_Vals_M0_fA_T1_chi2err.series',int2str(seriesnum*1000),'.slice',int2str(slice),'.mat'));
T10=[T10 Params.EstPrecontrast];
end
precontrastT1=median(T10)





% % .02 rest study
seriesnumStart=23
seriesnumEnd=25
% must also change filenames
scan_name='CV_Radial7Off_5ml('
MID=37

if doReconFlag==1
    load(['/v/raid1/' RawOwner '/MRIdata/Cardiac/' scanner '/' Patient '/RawData/meas_MID37_CV_Radial7Off_5ml_FID6341_Kspace.mat'],'-mat');
end

for seriesnum=seriesnumStart:seriesnumEnd
	slice=seriesnum-seriesnumStart+1
    template_dicom_directory=strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',Patient,'/DicomData/',scan_name,int2str(seriesnum),')')
    output_dicom_directory=strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',Patient,'/ReconData/Subsets/',scan_name );
    outname=strcat('MID_',int2str(MID),'_imgSW_slice_',int2str(slice),'_subset');  % will add on subset # and window size when use as filename and 

    if doReconFlag
        kSpaceOneSlice = single(kSpace(:,:,1:3:end,startFrame:end,slice) * 1e9);    
        % used as series description
        func_directFBPandShare_3subsets(kSpaceOneSlice, template_dicom_directory, output_dicom_directory, slice, seriesnum*1000, outname );
        % only one slice at a time for now!!!
    end
  
    if flipFlag
        for subset=1:3
            % get directory name
            dirname=strcat(output_dicom_directory,int2str(seriesnum*1000+subset),')/')
            files=dir(strcat(dirname,'*.dcm'))
            for bb=1:length(files) 
                tmpp=dicomread( strcat(dirname,files(bb).name));
                %tmpp=flipud(tmpp); 
                %imgrr(:,:,bb)=rot90(tmpp);
                imgrr(:,:,bb)=flipud(tmpp);
            end

            seriesdesc = outname;
            seriesdesc_out = strcat(seriesdesc,int2str(subset),'_window',int2str(10*modifyWindowSize));
            % write out as dicom:
            func_dicomHeader_OrderByTime_SkipFirst(template_dicom_directory, imgrr, [output_dicom_directory num2str(seriesnum*1000+subset) ')/' outname num2str(subset) '_window' int2str(10*modifyWindowSize) '_'], 1000*seriesnum+subset, seriesdesc_out);
        end
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
        fit_subsets(scanner, Patient, seriesnum*1000, slice, precontrastT1);
    end

end
   
%keyboard
   

% adeno study
seriesnumStart=28
seriesnumEnd=30
% must also change filenames
scan_name='CV_Radial7Off_adeno_7.5ml('
MID=40
if doReconFlag==1
    load(['/v/raid1/' RawOwner '/MRIdata/Cardiac/' scanner '/' Patient '/RawData/meas_MID40_CV_Radial7Off_adeno_7.5ml_FID6344_Kspace'],'-mat');
end


for seriesnum=seriesnumStart:seriesnumEnd
	slice=seriesnum-seriesnumStart+1
    template_dicom_directory=strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',Patient,'/DicomData/',scan_name,int2str(seriesnum),')')
    output_dicom_directory=strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',Patient,'/ReconData/Subsets/',scan_name );
    outname=strcat('MID_',int2str(MID),'_imgSW_slice_',int2str(slice),'_subset');  % will add on subset # and window size when use as filename and 

    if doReconFlag
        kSpaceOneSlice = single(kSpace(:,:,1:3:end,startFrame:end,slice) * 1e9);    
        % used as series description
        func_directFBPandShare_3subsets(kSpaceOneSlice, template_dicom_directory, output_dicom_directory, slice, seriesnum*1000, outname );
        % only one slice at a time for now!!!
    end
  
    if flipFlag
        for subset=1:3
            % get directory name
            dirname=strcat(output_dicom_directory,int2str(seriesnum*1000+subset),')/')
            files=dir(strcat(dirname,'*.dcm'))
            for bb=1:length(files) 
                tmpp=dicomread( strcat(dirname,files(bb).name));
                %tmpp=flipud(tmpp); 
                %imgrr(:,:,bb)=rot90(tmpp);
                imgrr(:,:,bb)=flipud(tmpp);
            end

            seriesdesc = outname;
            seriesdesc_out = strcat(seriesdesc,int2str(subset),'_window',int2str(10*modifyWindowSize));
            % write out as dicom:
            func_dicomHeader_OrderByTime_SkipFirst(template_dicom_directory, imgrr, [output_dicom_directory num2str(seriesnum*1000+subset) ')/' outname num2str(subset) '_window' int2str(10*modifyWindowSize) '_'], 1000*seriesnum+subset, seriesdesc_out);
        end
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
        fit_subsets(scanner,Patient, seriesnum*1000, slice, precontrastT1);
    end
end
   
   


% lexi study
seriesnumStart=55
seriesnumEnd=57
% must also change filenames
scan_name='CV_Radial7Off_lexi_7.5ml('
MID=56
if doReconFlag==1
    load(['/v/raid1/' RawOwner '/MRIdata/Cardiac/' scanner '/' Patient '/RawData/meas_MID56_CV_Radial7Off_lexi_7.5ml_FID6360_Kspace'],'-mat');
end


for seriesnum=seriesnumStart:seriesnumEnd
	slice=seriesnum-seriesnumStart+1
    template_dicom_directory=strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',Patient,'/DicomData/',scan_name,int2str(seriesnum),')')
    output_dicom_directory=strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',Patient,'/ReconData/Subsets/',scan_name );
    outname=strcat('MID_',int2str(MID),'_imgSW_slice_',int2str(slice),'_subset');  % will add on subset # and window size when use as filename and 

    if doReconFlag
        kSpaceOneSlice = single(kSpace(:,:,1:3:end,startFrame:end,slice) * 1e9);    
        % used as series description
        func_directFBPandShare_3subsets(kSpaceOneSlice, template_dicom_directory, output_dicom_directory, slice, seriesnum*1000, outname );
        % only one slice at a time for now!!!
    end
  
    if flipFlag
        for subset=1:3
            % get directory name
            dirname=strcat(output_dicom_directory,int2str(seriesnum*1000+subset),')/')
            files=dir(strcat(dirname,'*.dcm'))
            for bb=1:length(files) 
                tmpp=dicomread( strcat(dirname,files(bb).name));
                %tmpp=flipud(tmpp); 
                %imgrr(:,:,bb)=rot90(tmpp);
                imgrr(:,:,bb)=flipud(tmpp);
            end

            seriesdesc = outname;
            seriesdesc_out = strcat(seriesdesc,int2str(subset),'_window',int2str(10*modifyWindowSize));
            % write out as dicom:
            func_dicomHeader_OrderByTime_SkipFirst(template_dicom_directory, imgrr, [output_dicom_directory num2str(seriesnum*1000+subset) ')/' outname num2str(subset) '_window' int2str(10*modifyWindowSize) '_'], 1000*seriesnum+subset, seriesdesc_out);
        end
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
        fit_subsets(scanner, Patient, seriesnum*1000, slice, precontrastT1);
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



