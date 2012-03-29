% read in raw, recon with FBP and 3 subsets, and write out as dicom

% Note, views only shared within each 72 ray time frame (could extend to do 3 other time frames for 288 rays total). 
% Script needs dicom data unpacked by dicomsorter_recursive.m    
% 
% Takes in raw (.mat) k-space file and writes out recon with appropriate dicom headers, in 
% subdirectory to the Dicom data called Recon (can be altered. Note: don't write out to same Dicom directory,
% as subsequent uses of this program could get confused)

clear
scanner='Verio'
Patient='P052710_alt'
types ={'low','high','adeno','lexi'};
MIDs = [168    169    170     184];
for i=1:length(types)
    disp([types{i} ' : ' num2str(MIDs(i))]);
end
startFrame=1;
%addpath /Users/ed/Documents/ShareWithWindows/matlab/generalCodes
addpath /v/raid1/ed/src/matlab/Radial
addpath /v/raid1/ed/src/matlab/generalCodes
addpath /v/raid1/bmatthew/brian-Code
%addpath('/v/raid1/npack/Code/mpi2dCode_THK_092209/')

processingDir = ['/v/raid1/bmatthew/MRIdata/Cardiac/' scanner '/' Patient '/Processing'];
if(~exist(processingDir,'dir'))
    disp('That patient folder doesn''t exist, or the Processing directory hasn''t been created');
    return;
end
cd(processingDir);
cd('..')
if(exist('DicomData')~=7)
    disp('I don''t think you''re in the Processing folder or the dicom folder isn''t linked');
    return;
else
    dicomDataFolder = cd('DicomData');
    cd(processingDir);
    cd('..')
end
if(exist('ReconData') ~=7)
    disp('I don''t think you''re in the Processing folder or the Recon folder isn''t linked');
    return;
else
    ReconDataFolder = cd('ReconData');
    cd(processingDir);
    cd('..')
end
doReconFlag=1;  % in case recons already done, don't need to repeat. 
% what if mother not analyzed?
% do all of these - writes to my directory. 
% then cp -pr
% /v/raid1/ed/MRIdata/Cardiac/Verio/P041910/ReconData/*00[1-3]*
% /v/raid1/gadluru/MRIdata/Cardiac/Verio/P041910/ReconData/Subsets
doProcessing=1;  % only do this if "mother" series have been registered and segmented. 

% first injection, so we can get T1_0
cd(processingDir);
if(~exist('MIDinfolookup.mat','file'))
    [MIDinfolookup,SeriesNumberinfolookup] = findMIDrelations();
    if(isempty(MIDinfolookup))
        disp('I don''t think you''re in the Processing folder');
        return;
    end
else
    load('MIDinfolookup.mat');
end


for MIDi = 1:length(MIDs)
    MID = MIDs(MIDi);
    seriesStruct = MIDinfolookup(MID);

    seriesnumStart=min(seriesStruct.SeriesNumbers);
    seriesnumEnd=max(seriesStruct.SeriesNumbers);
    % must also change filenames
    scan_name=[seriesStruct.SeriesDescription '('];
    disp([num2str(MID) ' : ' scan_name num2str(seriesnumStart:seriesnumEnd) ')']);
    if doReconFlag==1
        cd('../RawData');
        matfileName = dir(['*MID' num2str(MID) '*Kspace*']);
        if(isempty(matfileName))
            disp(['Cannot find mat file with MID ' num2str(MID)]);
            return;
        end
        load(matfileName(1).name);
    end
    for seriesnum=seriesnumStart:seriesnumEnd
        slice=seriesnum-seriesnumStart+1;
        %if doReconFlag
        outname=strcat('MID_',int2str(MID),'_imgSW_slice_',int2str(slice),'_subset');  % will add on subset # and window size when use as filename and 
        % used as series description
        if(~exist([ReconDataFolder '/Subsets/' outname],'file'))
            kSpaceOneSlice = single(kSpace(:,:,1:3:end,startFrame:end,slice) * 1e9);
            %template_dicom_directory=strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',Patient,'/DicomData/',scan_name,int2str(seriesnum),')');
            template_dicom_directory=[dicomDataFolder '/' scan_name int2str(seriesnum) ')'];
            output_dicom_directory=[ReconDataFolder '/' scan_name];
            func_directFBPandShare_3subsetsED(kSpaceOneSlice, template_dicom_directory, output_dicom_directory, slice, seriesnum*1000, outname );
        end

        if doProcessing
            % get current directory so we can return to it afterwards
            cd(processingDir);
            for subset=1:3
                folder=strcat(scan_name,int2str(seriesnum*1000+subset),')');
                eval(['mpi2d stg=0.3 DoThisFolder=''' folder ''''])
                delete('mpi2d.par')
                if(~exist(['Output/blood_polyCoords.study' num2str(seriesnum*1000) '.slice' num2str(slice)],'file'))
                    disp(['Processing hasn''t been done on series ' num2str(seriesnum*1000)])
                    return;
                end
                parfile=strcat('rad.series',int2str(seriesnum*1000+subset),'.slice',int2str(slice),'.subset',int2str(subset),'.par');
                copyfile(parfile,'mpi2d.par')
                %copyfile('rad.series18002.slice1.subset2.par','mpi2d.par')
                mpi2d stg=[3.2 3.3]
            end
            if(MIDi == 1)
                precontrastT1=fit_subsets(Patient, seriesnum*1000, slice);
            else
                fit_subsets(Patient, seriesnum*1000, slice, precontrastT1);
            end
        end

    end
    if(MIDi==1)
        % to pull out precontrastT1 for other studies:
        T10=[];
        for seriesnum=seriesnumStart:seriesnumEnd
            slice=seriesnum-seriesnumStart+1;
            load(strcat(processingDir,'/Output/MultiSRT_Vals_M0_fA_T1_chi2err.series',int2str(seriesnum*1000),'.slice',int2str(slice),'.mat'));
            T10=[T10 Params.EstPrecontrast]; %#ok<AGROW>
        end
        precontrastT1=median(T10) %#ok<NOPTS>
    end
end








