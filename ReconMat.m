function ReconMat(MID, seriesNum,slice,OrientationTokens,coils)
                            %only specify JustaFew if you want a few time frames 
    disp('Collecting Coil opinions and adding Dicom Headers');

    %must be in the Processing folder
    Processingwd = pwd;
    %% set parameters
    showcoil=1;
    if(~exist('coils','var'))
        coils = [1 4 7 10];
    end
    
    
    %% load an example mat file just to see size
    cd('../ReconData/mat_files/combined/');
    regexFinderString = ['*MID' num2str(MID) '*.mat'];
    rawMatFiles = dir(regexFinderString);
    if(isempty(rawMatFiles))
        disp(['I couldn''t find any mat files for MID: ' num2str(MID)]);
        keyboard;
    end
    x = load(rawMatFiles(1).name);
    myfields = fields(x);
    if(length(myfields) > 1)
        disp('Mat files have more than 1 variable stored in them');
        return;
    end
    img_est = x.(myfields{1});
    [sx sy tmax]=size(img_est);
    if(tmax > 40)
        myreferenceFrame = 20;
    else
        myreferenceFrame = 6;
    end
    sos=zeros(size(img_est));
    % clear img_est
    if(~isempty(strfind(rawMatFiles(1).name,'Coil')) || coils~=1)
        %% seperate the mat files into slices
        mymatfiles = containers.Map({0},{''});remove(mymatfiles,0);
        for filei = 1:length(rawMatFiles)
            %is this the right slice number?
            if(~isempty(strfind(rawMatFiles(filei).name,['slice' num2str(slice)])))
                for coil=coils
                    %does it match one of the coils we have listed?
                    if(~isempty(strfind(rawMatFiles(filei).name,['Coil' num2str(coil) '_'])))
                        %Yay!!! we add it
                        mymatfiles(coil) = rawMatFiles(filei).name;
                        break;
                    end
                end
            end
        end
    else
        mymatfiles = containers.Map({0},{''});mymatfiles.remove(0);
        regexFinderString = ['*MID' num2str(MID) '*slice' num2str(slice) '*.mat'];
        rawMatFiles = dir(regexFinderString);
        if(length(rawMatFiles) > 1)
            disp(['There were more than 1 mat file that has that MID number(' num2str(MID) ') and slice number(' num2str(slice) ')']);
            disp('Did you mean to specify coil numbers, or did you forget to name the files with ''Coil#'' in the filename?');
            disp('If you meant to only tie 1 file to dicoms please don''t have Coil in the mat file names');
            return;
        end
        coils = 1;
        mymatfiles(1) = rawMatFiles(1).name;
    end
    
    if(length(coils) ~= 1)
        %% extract each coil opinion and merge with other coils and save out
        %% dicoms
        %get the raw file and accumulate it
        a = ceil(sqrt(length(coils)+1));  % move  12/12/10  EVRD
        b = ceil((length(coils)+1)/a);
        coilNumbers = coils;
        for coili = 1:length(coilNumbers)
            rawFile = mymatfiles(coilNumbers(coili));
            %load(rawFile);
            % need something here to check if loaded data is named img_est
            % EVRD  12/12/10, took from above
            x = load(rawFile);
            myfields = fields(x);
            img_est = x.(myfields{1});
            if(showcoil)
                figure(slice);
                subplot(a,b,coili+1), imagesc(img_est(:,:,myreferenceFrame)), colormap gray, axis image off;
                title(['unflipped Coil ' num2str(coilNumbers(coili))]);
                pause(.1);
            end

            %do the accumulation
            sos = sos + img_est.*conj(img_est);
        end
        sos = sqrt(sos);
    else
        for filei = 1:length(rawMatFiles)
            x = load(rawMatFiles(1).name);
            %load(rawFile);
            % need something here to check if loaded data is named img_est
            % EVRD  12/12/10, took from above
            myfields = fields(x);
            img_est = x.(myfields{1});
            if(showcoil)
                figure(slice);
                imagesc(img_est(:,:,myreferenceFrame)), colormap gray, axis image off;
                title(['unflipped slice ' num2str(filei)]);
                pause(.1);
            end
            sos = abs(img_est);
    end
    
    
    
    imgrr = sos((sx/4+1):(sx-sx/4),(sy/4+1):(sy-sy/4),:);

    %% flippings
    if(exist('OrientationTokens','var') && ~strcmp(OrientationTokens{1},''))
        clear newimgrr
        for t=1:size(imgrr,3)
            for opperationi = 1:length(OrientationTokens)
                operation = OrientationTokens{opperationi};
                if(strcmp(operation,'flipud'))
                    newimgrr(:,:,t) = flipud(imgrr(:,:,t)); %#ok<*AGROW>
                end
                if(strcmp(operation,'fliplr'))
                    newimgrr(:,:,t) = fliplr(imgrr(:,:,t));
                end
                if(strcmp(operation,'rot90'))
                    newimgrr(:,:,t) = rot90(imgrr(:,:,t));
                end
                if(strcmp(operation,'rot-90'))
                    newimgrr(:,:,t) = rot90(imgrr(:,:,t),-1);
                end
            end
        end
        imgrr = newimgrr;
    end
    flippings = '';
    if(exist('OrientationTokens','var'))
        for i=1:length(OrientationTokens)
            flippings = [flippings ' - ' OrientationTokens{i}];
        end
    end
    a = ceil(sqrt(length(coils)+1));  % move  12/12/10  EVRD
    b = ceil((length(coils)+1)/a);
    subplot(a,b,1), imagesc(imgrr(:,:,myreferenceFrame)), colormap gray, axis image off;
    title(['Kosher Orientation Slice: ' num2str(slice)  ' { ' flippings ' }']);
    pause(.1);

    %% get the dicom folder and put the dicom headers on these images
    cd(Processingwd); cd('../DicomData');
    targetReconDicomFolder = dir(['*(' num2str(seriesNum) ')']);
    targetReconDicomFolder = targetReconDicomFolder(1).name;
    targetDicomDataFolder = targetReconDicomFolder;
    targetReconDicomFolder = strrep(targetReconDicomFolder,num2str(seriesNum) ,num2str(seriesNum*1000));
    cd(Processingwd);  %EVRD
    cd('../ReconData');
    %delete previous dicoms
    if(exist(targetReconDicomFolder,'dir') == 7)
        cd(targetReconDicomFolder);
        delete('*.dcm');        %he's just sampling a single time frame
    else
        mkdir(targetReconDicomFolder);
        cd(targetReconDicomFolder);
    end
    ReconDatawd = pwd;   %this is where to place the reconstructed dicoms
    cd([Processingwd '/../DicomData/' targetDicomDataFolder]);
    dicomFolder = pwd;  %this is the absolute path to the scanner dicom directory
    dicomFolder = [dicomFolder '/'];
    func_dicomHeader_OrderByTime_SkipFirst(dicomFolder, imgrr, [ReconDatawd '/'], seriesNum*1000, strrep(targetDicomDataFolder,['(' num2str(seriesNum) ')'],''));
    cd(Processingwd);
end