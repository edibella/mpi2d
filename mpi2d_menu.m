function mpi2d_menu()
global menuv
temp = pwd;
menuv.pwd = temp;
if(strcmp(temp((end-9):end),'Processing') || strcmp(temp((end-8):end),'DicomData') || strcmp(temp((end-8):end),'ReconData'))
    cd('..');
end

%which folder should I be looking in
if(isdir('ReconData'))
    menuv.dicomfolder = 'ReconData';
elseif(isdir('DicomData'))
    menuv.dicomfolder = 'DicomData';
else
    disp('Could not find the DicomData/ReconData folder');
    return;
end

%what files have already been processed
cd('Processing');
    parfiles = dir('*.par');
cd('..');

cd(menuv.dicomfolder);
    movieFoldercandidates = dir('*');
    menuv.goodMovieFolders = {};
    clear temp;
    temp.hasParfile = 0;
    temp.hascinefile = 0;
    temp.hascontourfile = 0;
    temp.hasbeenregistered = 0;
    temp.hasKflowfile = 0;
    menuv.Movies = containers.Map({movieFoldercandidates(1).name},{zeros(2,2)});
    menuv.parfiles = containers.Map({movieFoldercandidates(1).name},{'dummy'});
    menuv.identifier = containers.Map({movieFoldercandidates(1).name},{[0 0]});
    menuv.fileExistances = containers.Map({movieFoldercandidates(1).name},{temp});
    menuv.goodMoviethumbnail = containers.Map({movieFoldercandidates(1).name},{zeros(2,2)});   % this first one is guaranteed to be . but is of the correct type for our map
    for i=3:length(movieFoldercandidates)

        %skip all files
        if(~movieFoldercandidates(i).isdir)
            continue;
        end
        cd(char(movieFoldercandidates(i).name));
            possibleDicomFiles = dir('*.dcm');
            %skip all folders that are not full length movies
            if(length(possibleDicomFiles) < 30)
                cd('..');
                continue;
            end
            disp(['Found: ' char(movieFoldercandidates(i).name)]);

            %pick 5 random frames and average them together as a thumbnail
            selection = ceil(rand(1,5)*length(possibleDicomFiles));
            clear thisThumbnail;
            for j=1:length(selection)
                dcmfile = char(possibleDicomFiles(selection(j)).name);
                dicomimage = dicomread(dcmfile);
                thisThumbnail(j,:,:) = dicomimage;
            end

            header = dicominfo(dcmfile);  %we'll just use the last valid dicom file
            %get the slice number
            PathName = char(movieFoldercandidates(i).name);
            beginningParen = strfind(PathName,'(');
            pattern = PathName(1:beginningParen(end));
            cd('..');
            matches = dir([pattern '*']);
            cd(char(movieFoldercandidates(i).name));
            tlist = cell([length(matches),1]);
            for ii=1:length(matches)
                tlist(ii) = {matches(ii).name};
            end
            matches = tlist;
            matches = sort(matches);
            slice = 1;
            for ii=1:size(matches)
                index = strfind(PathName,'/');
                if(strcmp(char(matches(ii)),PathName) > 0)
                    slice = ii;
                    break;
                end
            end
            disp(['I think I am slice# ' num2str(slice)]);

            %save away the numbers used to construct the filenames of stuff
            %needed
            menuv.identifier(movieFoldercandidates(i).name) = [slice header.SeriesNumber];
            cd('..'),cd('..'), cd('Processing');
            clear temp;
            %construct the parfilename and check to see if it exists
            temp.parfilename = ['series' int2str(header.SeriesNumber) '.slice' int2str(slice) '.par'];
            temp.hasParfile = exist( temp.parfilename)>0;
            cd('Output');
                temp.cinefilename = ['cinemri1.study' int2str(header.SeriesNumber) '.slice',int2str(slice) '.mat'];
                temp.contourfilename = ['blood_polyCoords.study' int2str(header.SeriesNumber) '.slice',int2str(slice)];
                temp.registeredfilename = ['shiftsMAN.study' int2str(header.SeriesNumber) '.slice',int2str(slice) '.txt'];
                temp.Kflowfiles = dir(['flowvalues.study' int2str(header.SeriesNumber) '.slice' int2str(slice) '*']);
                temp.hascinefile = exist(temp.cinefilename)>0;
                temp.hascontourfile = exist(temp.contourfilename)>0;
                temp.hasKflowfile = length(temp.Kflowfiles)> 0;
                temp.hasbeenregistered = exist(temp.registeredfilename)>0;
                
                %if the cine movie exists, use that for the snapshot
                if(temp.hascinefile)
                    load(temp.cinefilename);
                    clear thisThumbnail;
                    for j=1:5
                        thisThumbnail(j,:,:) = cinemri1(:,:,ceil(rand(1,1)*size(cinemri1,3)));
                    end
                end
            cd('..');
            menuv.fileExistances(movieFoldercandidates(i).name) = temp;
        cd('..'), cd(menuv.dicomfolder);
        
        %save the thumbnail away in a map
        menuv.goodMoviethumbnail(movieFoldercandidates(i).name) = squeeze(mean(thisThumbnail,1));
        menuv.goodMovieFolders = vertcat(menuv.goodMovieFolders,movieFoldercandidates(i).name);
    end
cd('..');
menuv.s = ceil(sqrt(length(menuv.goodMovieFolders)));
menuv.subplotToFolderName = containers.Map({0},{menuv.goodMovieFolders{1}});
h = figure(67);
ontop = 1;
for i=1:length(menuv.goodMovieFolders)
    menuv.h(i) = subplot(menuv.s,menuv.s,i);
    menuv.subplotToFolderName(menuv.h(i)) = menuv.goodMovieFolders{i};
    imagesc(menuv.goodMoviethumbnail(menuv.goodMovieFolders{i})), colormap gray,axis off;
    temp = menuv.fileExistances(menuv.goodMovieFolders{i});
    parexists = ['P:' num2str(temp.hasParfile) ' '];
    Registered = ['R:' num2str(temp.hasbeenregistered) ' '];
    contoursExist = ['C:' num2str(temp.hascontourfile) ' '];
    kFlowExist = ['K:' num2str(temp.hasKflowfile) ' '];
    ontop = bitand((mod(i,menuv.s+1)) , 1);
    if(ontop) 
       list = {strrep(menuv.goodMovieFolders{i},'_',''),' ',[ parexists Registered contoursExist kFlowExist]};
    else
        list = {' ',strrep(menuv.goodMovieFolders{i},'_',''),[ parexists Registered contoursExist kFlowExist]};
    end
    
    title(list);
end
set(h,'WindowButtonDownFcn',{@f_wbfcn});
end




function [] = f_wbfcn(hObject, eventdata, handles)
global menuv
%right or left click
mouseside=get(gcf,'SelectionType');

%get the subplot handle they clicked on
if(~isKey(menuv.subplotToFolderName,gca))
    return;
end
folderName = menuv.subplotToFolderName(gca);

%find what index so we can do subplot(x,x,index)
for i=1:length(menuv.goodMovieFolders)
    if(strcmp(menuv.goodMovieFolders{i},folderName))
        subplotindex = i;
        break;
    end
end
fileExistances = menuv.fileExistances(menuv.goodMovieFolders{subplotindex});

if(strcmp(mouseside,'alt'))    % right click
    disp('------------------------------------------');
    if(fileExistances.hasParfile) disp(['Parfile= ' fileExistances.parfilename]); end
    if(fileExistances.hascinefile) disp(['CineFilename= ' fileExistances.cinefilename]); end
    if(fileExistances.hasbeenregistered) disp(['Shift filename= ' fileExistances.registeredfilename]); end
    if(fileExistances.hasKflowfile) 
        disp(['K-flow values: ']); 
        cd([menuv.pwd '/Processing/Output']);
        for i=1:length(fileExistances.Kflowfiles)
            disp(['   ' char(fileExistances.Kflowfiles(i).name)]);
            fid = fopen(char(fileExistances.Kflowfiles(i).name));
            tline = fgetl(fid);
            while ischar(tline)
                if(~isempty(strfind(tline,'Ktrans')))
                    disp(['      ' tline]);
                end
                if(~isempty(strfind(tline,'flowMean')))
                    disp(['      ' tline]);
                end
                tline = fgetl(fid);
            end

            fclose(fid);
        end
        cd('../..');
    end
    disp(['1: Process(Auto-Register + par file creation)']);
    disp(['2: create/Fix contours + Manual registration tweaks']);
    disp(['3: Process contours and get k-values curves']);
    disp(['4: Unprocess all associated files']);
    key = input('');
    switch key
        case 1
            cd([menuv.pwd '/Processing']);
            temp = ['mpi2d stg=0.2 DoThisFolder=' folderName];
            eval(temp);
        case 2
            cd([menuv.pwd '/Processing']);
            copyfile(fileExistances.parfilename,'mpi2d.par');
            mpi2d stg=3.11
            fileExistances.contoursExist = 1;
        case 3
            cd([menuv.pwd '/Processing']);
            copyfile(fileExistances.parfilename,'mpi2d.par');
            mpi2d stg=[3.2 3.3 3.4 4.1 4.2]
            fileExistances.kFlowExist = 1;
        case 4
            cd([menuv.pwd '/Processing']);
            A = menuv.identifier(folderName);
            slice = A(1); series = A(2);
            if(fileExistances.hasParfile) 
                %delete( fileExistances.parfilename); 
            end
            filesToDelete = dir(['*.study' int2str(series) '.slice' int2str(slice) '*']);
            for i=1:length(filesToDelete)
                delete( char(filesToDelete(i).name));
            end
            cd('Output');
                filesToDelete = dir(['*.study' int2str(series) '.slice' int2str(slice) '*']);
                for i=1:length(filesToDelete)
                    delete(char(filesToDelete(i).name));
                end
                fileExistances.hasParfile = 0;
                fileExistances.hascinefile = 0;
                fileExistances.hasbeenregistered = 0;
                fileExistances.hasKflowfile = 0;
                fileExistances.hascontourfile = 0;
            cd('..');
        otherwise
            disp('?');
    end
    figure(67)
    subplot(menuv.s,menuv.s,subplotindex);
    parexists = ['P:' num2str(fileExistances.hasParfile) ' '];
    Registered = ['R:' num2str(fileExistances.hasbeenregistered) ' '];
    contoursExist = ['C:' num2str(fileExistances.hascontourfile) ' '];
    kFlowExist = ['K:' num2str(fileExistances.hasKflowfile) ' '];
    ontop = bitand((mod(subplotindex,menuv.s+1)) , 1);
    if(ontop) 
       list = {strrep(folderName,'_',''),' ',[ parexists Registered contoursExist kFlowExist]};
    else
        list = {' ',strrep(folderName,'_',''),[ parexists Registered contoursExist kFlowExist]};
    end
    title(list);
elseif(strcmp(mouseside,'normal'))   %left click
    %check to see if we've pulled the dicom files out before
    if(~isKey(menuv.Movies,folderName))
        %check to see if the cinemri file exists
        if(fileExistances.hascinefile)
            cd([menuv.pwd '/Processing/Output']);
            load(fileExistances.cinefilename);
            cinemri = cinemri1;
            cd('../..');
        else

            %else load it up from the disk
            h = waitbar(0,'Reading Dicomfiles');
            cd([menuv.pwd '/' menuv.dicomfolder '/' folderName]);
                dcmfiles = dir('*.dcm');
                sample = dicomread(dcmfiles(1).name);
                cinemri = zeros(size(sample,1),size(sample,2),length(dcmfiles));
                for i=1:length(dcmfiles)
                    waitbar(i/length(dcmfiles));
                    header = dicominfo(dcmfiles(i).name);
                    index = header.InstanceNumber;
                    dcmimage = dicomread(dcmfiles(i).name);
                    cinemri(:,:,index) = dcmimage;
                end
            cd('..');cd('..');
            close(h);
        end
        menuv.Movies(folderName) = cinemri;
    end

    %play the movie
    cinemri = menuv.Movies(folderName);
    for i=1:size(cinemri,3)
        subplot(menuv.s,menuv.s,subplotindex);
        imagesc(squeeze(cinemri(:,:,i))), colormap gray,axis off, pause(.05);
    end
    parexists = ['P:' num2str(fileExistances.hasParfile) ' '];
    Registered = ['R:' num2str(fileExistances.hasbeenregistered) ' '];
    contoursExist = ['C:' num2str(fileExistances.hascontourfile) ' '];
    kFlowExist = ['K:' num2str(fileExistances.hasKflowfile) ' '];
    ontop = bitand((mod(subplotindex,menuv.s+1)) , 1);
    if(ontop) 
       list = {strrep(folderName,'_',''),' ',[ parexists Registered contoursExist kFlowExist]};
    else
        list = {' ',strrep(folderName,'_',''),[ parexists Registered contoursExist kFlowExist]};
    end
    title(list);
end
menuv.fileExistances(menuv.goodMovieFolders{subplotindex}) = fileExistances;
end