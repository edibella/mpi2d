%
% Version 6/2004
%
% mpi2d lst=deliver
% produces mpi2d.m in the current directory

% 
% to convert to *.p code use
% pcode mpi2d.m
%
function maincode(varargin) %#ok<FNDEF>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   CODE HEADER: HELP MESSAGES AND INPUT OF PARAMETERS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
version = 1.5; %#ok<NASGU>

tic;   % starts stopwatch timer, call toc to get time 
code='mpi2d';     % myocardial perfusion imaging
data=[code '.mat'];

if(nargin==0) % help message printout mode only

 disp([ ' ' code ' : ']);
 disp('Myocardal Perfusion Imaging 2D (mpi2d) for dynamic contrast MRI');
 dispHelp(code);
 disp('Contact:'); disp('Ed DiBella ed@ucair.med.utah.edu');
 disp(['Version: ' codedate ]);
  return;

else          % get parameters, do the work
 p={' ',' '};
end

p  = buf2par(p,arg2buf(nargin,varargin));
lst = getpar(p,'char','lst',   '-1','List the names of stages');

switch lower(lst)
 case 'par'
  p=[];  
 case 'deliver'
  disp('Producing delivery code');
  disp('Disabled, make backup and check before using');
  deliver(code,'../code');
  disp(['Version: ', codedate ]);
 case 'all'
  lst=[0:.001:10]; %#ok<NBRAK>
 otherwise
  lst=txt2double(lst);
end

in = getpar(p,'char',  'in',[ code '.par'],'Input parameter file name');
p  = buf2par(p,readtxt(in));
stg = getpar(p,'double','stg',   '1','List of stages to execute');
wordy = getpar(p,'double','wordy','4'   ,'Messages: 0-silience, 4-very wordy');

% Get parameters from file

infile = getpar(p,'char','infile','cine.dat','Raw data input file');
coilsensfile = getpar(p,'char','coilsensfile','','Raw data coil sens file'); %#ok<NASGU>
inendi = getpar(p,'char','inendi','ieee-be','Input file endianess');
intype = getpar(p,'char','intype','ushort','Input file word type');
inpnxy = getpar(p,'double','inpnxy','[256 256]','Input image size x, y');
rangex = getpar(p,'double','rangex','50:180','Input window range in x');
HorzOffset = getpar(p,'double','HorzOffset','0','Offset for text in stage 6.1');
Alpha = getpar(p,'double','Alpha','.5','Alpha channel for k-trans overlay in stage 6.1');
MyoFatness = getpar(p,'double','MyoFatness','1','How fat the blending should be for the k-trans overlay in stage 6.1');
FinalShading = getpar(p,'double','FinalShading','1','In stage 6.1 if the final display should be shaded');
flipx = getpar(p,'double','flipx','0','Flip x direction');
flipy = getpar(p,'double','flipy','0','Flip y direction');
myothreshold = getpar(p,'double','myothreshold','.3','threshold for the myocardium');
rotNinety = getpar(p,'double','rotNinety','0','rotate 90');
rangey = getpar(p,'double','rangey','50:180','Input window range in y');
ranget = getpar(p,'double','ranget','1:20','Input window range in t');
framesToSelect = getpar(p,'double','framesToSelect','[]','Frames to leave out due to cardiac phase, NOT REPLACED');
framesToSkip = getpar(p,'double','framesToSkip','[]','Frames to leave out due to motion, will be replaced by linear interp of neighbors');
delta_t = getpar(p,'double','delta_t','2','Time between frames (same slice), in seconds');
timeStampFile = getpar(p,'char','timeStampFile','','Name and path of .mat file with timestamps');
timeStampFileAIF = getpar(p,'char','timeStampFileAIF','','Name and path of .mat file with timestamps for low dose AIF');
flagTimeStamps=0; %#ok<NASGU>
if exist(timeStampFile)
	flagTimeStamps=1; %#ok<NASGU>
end
sliceNum = getpar(p,'double','sliceNum','1','Slice number being processed');
seriesNum = getpar(p,'double','studyNum','-1','Series number being processed');
if seriesNum==-1
   seriesNum = getpar(p,'double','seriesNum','-1','Series number being processed');
end
% don't want to use default, leaving studyNum for historical compability

numPreContrast_bld = getpar(p,'double','numPreContrast_bld','5','self-explanoat');
numPreContrast_tiss = getpar(p,'double','numPreContrast_tiss','6','self-explanoat');
numSkip = getpar(p,'double','numSkip','0','num frames to skip when converting to Gd, if intial glitches present');
numSkipUpslope = getpar(p,'double','numSkipUpslope','2','this'); %#ok<NASGU> %so don't compute upslope of some initial small bump  should code way of doig this
numSkipRegistration = getpar(p,'double','numSkipRegistration','0','num frames to skip when doing the registration, so can do initial frames manually'); %#ok<NASGU>

diffthld= getpar(p,'double','diffthld','1','Mask treshold, in std units'); %#ok<NASGU>
Niter= getpar(p,'double','Niter','1','Number of iterations for correction'); %#ok<NASGU>

delay = getpar(p,'double','delay','0.25','Frame delay in sec. for movie, stg 1.2'); %#ok<NASGU>

outfile = getpar(p,'char','outfile','cinemri.bin','Raw data output file'); %#ok<NASGU>
%staticonly = getpar(p,'double','staticonly','0','Static correction only=1 yes 0 no');  %#ok<NASGU> % meaning shifts only, no warping

numAzimuthalRegions = getpar(p,'double','numAzimuthalRegions','8','Choose number of regions, total number of regions will double if endoepiflag=1');
numRadialRegions = getpar(p,'double','numRadialRegions','1','Choose number of radial regions, 2 will be like endo and epi');
flagPixelwise = getpar(p,'double','flagPixelwise','0','Choose ROIs or pixelwise fits between contours');
% SHould add checks here that is only 0 or 1...
flagUseModel = getpar(p,'double','flagUseModel','0','Choose to use init_image rather than calculate it for deltaSI. Not yet supported for [Gd]. Added 4/04'); %#ok<NASGU>

flagUpsampleImage= getpar(p,'double','flagUpsampleImage','0','Upsample image by 2x in each direction with cubic interpolation. Or use 3.01 to do this after registration, for time reasons ');

referenceFrame = getpar(p,'double','referenceFrame','4','Frame to register all other frames to in stg 2.1');
useDeltaSI = getpar(p,'double','useDeltaSI','1','Instead of estimating gad. concentration, use signal differences');
useIntegralLinearFit = getpar(p,'double','useIntegralLinearFit','0','Instead of standard nonlinear fmincon');

pixelsizeX = getpar(p,'double','pixelsizeX','1','For accurate azimuthal angle division if pixels are not isotropic i.e. same size in x and y');
pixelsizeY = getpar(p,'double','pixelsizeY','1','For accurate azimuthal angle division if pixels are not isotropic i.e. same size in x and y');
pixelSizeRatio=pixelsizeX/pixelsizeY; %#ok<NASGU>
%     sliceNumAIF=5;  seriesNumAIF=14;
fixedVp = getpar(p,'double','fixedVp','99','Will estimate Vp for each region by default. Or will use given value for all regions.'); %#ok<NASGU>
fixedDelay = getpar(p,'double','fixedDelay','99','Will estimate time delay for each region by default. Or will use given time delay for all regions.');
seriesNumAIF = getpar(p,'double','seriesNumAIF','-1','Series number for AIF');
sliceNumAIF = getpar(p,'double','sliceNumAIF','-1','Slice number foor AIF');
scaleAIF = getpar(p,'double','scaleAIF','-1','Scale for AIF');
lastFrame = getpar(p,'double','lastFrame','0','last Frame for truncated models');
copyFromSeriesNum = getpar(p,'double','copyFromSeriesNum','0','study to copy from'); %#ok<NASGU>
DoThisFolder = getpar(p,'char','DoThisFolder','','which folder to process');
if(isempty(p)) % parameter printout mode only
  return;
else
  clear p
  if(lst(1)>=0), stg=[wordy 0 lst];
  else          stg=[wordy 1 stg]; end;
end
try
    if(exist('mpi2d.log'))
        fid = fopen('mpi2d.log','a'); 
    else
        fid = fopen('mpi2d.log','w'); 
    end
    fprintf(fid,'%s',['mpi2d(ver ' num2str(version) ') ran at ' num2str(date) ' with the following params(' num2str(varargin{:}) ') on seriesNum: ' num2str(seriesNum) ' and sliceNum: ' num2str(sliceNum) '/n']); 
    fclose(fid);
catch ME %#ok<NASGU>
    
end

if ~exist('Output','dir')
  [SUCCESS,MESSAGE,MESSAGEID]= mkdir('Output'); %#ok<ASGLU,NASGU>  % it claims it fails, but seems to work!
end
outpath='Output/';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   STAGES BEGIN
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ran(stg,1.1,'Overhead creates properly named dicom files');
    if(exist('Processing') ~=7)
        mkdir('Processing');
    end
    if(exist('DicomData') ~= 7)
        mkdir('DicomData');
    end
    cd('DicomData');
    disp('Please navigate to the top folder which contains the unsorted dicom files');
    if(strcmp(DoThisFolder,''))
        [PathName] = uigetdir('Please navigate to the top dicom folder');
        if(strcmp(PathName,''))
            return;
        end
    else
        PathName = DoThisFolder;
    end
    [outname] = pwd;
    
    %dicomsorter_recursive(PathName,outname);	 % to convert the Dicom Images to .mat files
    sortDicoms(PathName,outname);	 % to convert the Dicom Images to .mat files
    cd('..');
    if(exist('Output') == 7)
        cd('Output');
        files = dir('*');
        if(length(files) <3)
            cd('..');
            rmdir('Output');
        end
    end
    cd('Processing');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;if ran(stg,1.2,'Preps Folder to process and creates par files in preperation for processing');
    if(exist('Processing'))
        temp = pwd;
        if(isempty(strfind(temp,'Processing')))
            cd('Processing');
        end
    else
        if(exist('../Processing'))
            cd('../Processing')
        else
            if(exist('../../Processing'))
                cd('../../Processing');
            else
                disp('I don''t know where you are.  Please go to either the patient folder or preferably the Processing folder of the Patient');
                return;
            end
        end
    end
    myprocessingdir = pwd;
    if(exist('Output') ~= 7)
        mkdir('Output');
    end
    if(exist('../ReconData') == 7)
        cd('../ReconData');
    else
        cd('../DicomData');
    end
    
    dicompwd = pwd;
    candidateFolders = dir('*');
    if(strcmp(DoThisFolder,''))
        FoldersToDo = {};
        for i=1:length(candidateFolders)
            if(~candidateFolders(i).isdir), continue; end
            dicomFiles = dir([candidateFolders(i).name '/*.dcm']);
            if(length(dicomFiles) < 39)
                disp([candidateFolders(i).name 'has less than 39 dicoms.  Skipping']);
                continue;
            end

            FoldersToDo = vertcat(FoldersToDo,candidateFolders(i).name);
        end
    else
        FoldersToDo = {DoThisFolder};
    end
    
    dicomfiles = dir([FoldersToDo{1} '/*.dcm']);
    temp = dicomread([FoldersToDo{1} '/' dicomfiles(1).name]);
    wasEmpty = isempty(temp);
    if(wasEmpty)
        FoldersToDo = FoldersToDo(2:end);
    end
    index = 1;
    while(wasEmpty)
        if(index > length(FoldersToDo))
            disp('None of the dicom folders had readable dicoms.   aborting');
            return;
        end
        dicomfiles = dir([FoldersToDo{index} '/*.dcm']);
        temp = dicomread([FoldersToDo{index} '/' dicomfiles(1).name]);
        wasEmpty = isempty(temp);
        if(wasEmpty)
            FoldersToDo = FoldersToDo(2:end);
        end
    end
    
    %process all folders in FoldersToDo
    
    %% do the processing
    clear data 
    data.value = 0;
    %runningRV = []; runningLV = [];
    
    runningRV = containers.Map({''},{data}); remove(runningRV,'');
    runningLV = containers.Map({''},{data}); remove(runningLV,'');
    centerRV = containers.Map({''},{[]}); remove(centerRV,'');
    centerLV = containers.Map({''},{[]}); remove(centerLV,'');
    usePreviousRVLV = 0;
    MyVal = 0;   %really convoluted way of redoing the processing if there are a few slices that are unacceptable.  The redo is centered with the median LVRV position
    mymap = containers.Map({''},{data});
    remove(mymap,'');
    
    askedForParFileOverwrite = 0;
    parFileMap = containers.Map({''},{''}); parFileMap.remove('');
    indexOfFoldersToDo = 1:length(FoldersToDo);
    figure(1), clf;
    while(MyVal >= 0)
        a = ceil(sqrt(length(FoldersToDo)));
        b = ceil(length(FoldersToDo)/a);
        for i=indexOfFoldersToDo
            cd(dicompwd)
            cd(FoldersToDo{i});
            disp('--------------------');
            disp(num2str(i));
            disp(FoldersToDo{i});
            clear data
            dicoms = dir('*.dcm');
            [sx sy] = size(dicomread(dicoms(1).name));
            st = length(dicoms);
            %% read in dicoms and upsample
            nonUpsampled = zeros(sx,sy,st);
            time = zeros(st,1);
            header = dicominfo(dicoms(ceil(st/2)).name);
            infile = pwd;
            subplot(a,b,i),imagesc(dicomread(dicoms(1).name)), colormap gray, 
            if(mod(i,2) == 0)
                title({[strrep(header.SeriesDescription,'_','\_') ]}); %#ok<NBRAK>
            else
                title({[strrep(header.SeriesDescription,'_','\_')],''}); %#ok<NBRAK>
            end
            pause(.1);
            for j=1:length(dicoms)
                nonUpsampled(:,:,j) = dicomread(dicoms(j).name);
                header = dicominfo(dicoms(j).name);
                timestr = header.AcquisitionTime;
                time(j) = str2num(timestr(1:2))*60*60 + str2num(timestr(3:4))*60 + str2num(timestr(5:end));
                seriesNum = header.SeriesNumber;
            end
            [time IX] = sort(time);
            time = time-time(1);
            nonUpsampled = nonUpsampled(:,:,IX);
            %find the slice number
            cd(dicompwd);
            temp = FoldersToDo{i}(1:strfind(FoldersToDo{i},'('));
            injectionName = temp;
            base = FoldersToDo{i}; 
            if(~isempty(strfind(base,'lice')))
                index = strfind(base,'lice');
                base(index:(index+4)) = [];
                base(index-1) = '*';
            end
            results = dir([temp '*']);
            sliceNum = 0;
            for j=1:length(results)
                sliceNum = sliceNum + 1;
                if(strcmp(results(j).name,base))
                    break
                end
            end
            %% upsample
            temp = mpi_upsampleImages(nonUpsampled(:,:,1));
            [sx sy] = size(temp);
            fullcinemri1 = zeros(sx,sy,st);
            for t=1:st
                fullcinemri1(:,:,t) = mpi_upsampleImages(nonUpsampled(:,:,t));
            end
            %% pick the LVRV
            if(~usePreviousRVLV)
                [RV,LV] = FindLVRV(fullcinemri1,0);
            else
                RV = centerRV(injectionName);
                LV = centerLV(injectionName);
            end
            if(~isKey(runningRV,injectionName))
                clear temp
                temp.point = [];
                temp.index = [];
                runningRV(injectionName) = temp;
                runningLV(injectionName) = temp;
            end
            temp = runningRV(injectionName);temp.point = [temp.point;RV];temp.index = [temp.index i]; runningRV(injectionName) = temp;
            temp = runningLV(injectionName);temp.point = [temp.point;LV];temp.index = [temp.index i]; runningLV(injectionName) = temp; 
            windowWidth = sx/2.5;
            windowHeight = sy/2.5;
            rangex = round(LV(1) - windowWidth/2):round(LV(1) + windowWidth/2);
            rangey = round(LV(2) - windowHeight/2):round(LV(2) + windowHeight/2);
            RV(1) = max(RV(1),min(rangex)+1);
            RV(2) = max(RV(2),min(rangey)+1);
            temp = zeros(size(fullcinemri1,1),size(fullcinemri1,2),3);
            try
                temp2 = fullcinemri1(:,:,22);  
            catch
                temp2 = fullcinemri1(:,:,3);  % added 4/27/11 evrd
            end
            
            temp2 = temp2-min(temp2(:));
            temp2 = temp2/max(temp2(:));
            temp(:,:,1) = temp2;
            temp(:,:,2) = temp2;
            temp(:,:,3) = temp2;
            subplot(a,b,i), imagesc(temp), colormap gray, hold on;
            if(mod(i,2) == 0)
                title({[strrep(header.SeriesDescription,'_','\_') ]}); %#ok<NBRAK>
            else
                title({[strrep(header.SeriesDescription,'_','\_')],''}); %#ok<NBRAK>
            end
            plot(LV(2), LV(1),'ro',RV(2), RV(1),'o');
            text(LV(2)+10, LV(1)-15,'LV','Color',[1 1 1]);
            text(RV(2)+10, RV(1)-15,'RV','Color',[1 1 1]);
            plot([min(rangey) max(rangey) max(rangey) min(rangey) min(rangey)],[min(rangex) min(rangex) max(rangex) max(rangex) min(rangex)]);
            RV = RV - [min(rangex) min(rangey)]; %#ok<NASGU>
            LV = LV - [min(rangex) min(rangey)]; %#ok<NASGU>
            cd(myprocessingdir);   %go back
            %create the par file
            ParFileName = ['series' num2str(seriesNum) '.slice' num2str(sliceNum) '.par'];
            parFileMap(ParFileName) = header.SeriesDescription;
            if(exist(ParFileName) && ~askedForParFileOverwrite)
                choice = questdlg('I see some par files already here.  Should I use the rangex/y in them or overwrite with auto-cropping?','Par file usage','Reuse','Overwrite','Reuse');
                askedForParFileOverwrite = 1;
                switch choice
                    case 'Reuse'
                        ReadPar;
                    case 'Overwrite'
                end
            end
            if(~exist('Template.par'))
                disp(['I cannot find the Template.par file.  Please put it in a place I can find like (' myprocessingdir ') if you don''t have write permissions to the Code directory']);
                keyboard;
            end
            copyfile(which('Template.par'),ParFileName);
            clear temp;  
            temp.rangex = [num2str(min(rangex)) ':' num2str(max(rangex))];
            temp.rangey = [num2str(min(rangey)) ':' num2str(max(rangey))];
            temp.infile = [infile '/'];
            temp.studyNum = seriesNum;
            temp.sliceNum = sliceNum;
            temp.seriesNumAIF = seriesNum;
            temp.sliceNumAIF = sliceNum;
            temp.ranget = ['7:' num2str(st)];
            temp.lastFrame = st;
            temp.framesToSelect=['1:' num2str(st)];  % added 9/20/11, need
            %to check:  
            temp.timeStampFile = ['timeStampSer' num2str(seriesNum) '.mat'];
            updateParFile(ParFileName,temp);

            %create the timestamp file
            timeStamp = time;
            save(temp.timeStampFile,'timeStamp');

            
            save(['Output/RVLVpos.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.mat'],'RV','LV');
            if(exist('../ReconData') == 7)
                cd('../ReconData');
            else
                cd('../DicomData');
            end
            %mymap(FoldersToDo{i}) = data;
        end
        cd(myprocessingdir);
        if(MyVal == 1)
            MyVal = -1;
        end
        if(MyVal == 0)
            %sanity check that the RV and LV are about in the same place
            mykeys = keys(runningLV);
            flagged = [];
            for j=1:length(mykeys)
                injectionName = mykeys{j};
                centerLV(injectionName) = round(median(runningLV(injectionName).point,1));
                clear dist
                for i=1:size(runningLV(injectionName).point,1)
                    dist(i) = norm(centerLV(injectionName) - runningLV(injectionName).point(i,:));
                end
                meandist = mean(dist);
                %stddist = std(dist);
                stddist = 20;
                flaggedLV = runningLV(injectionName).index(dist>(stddist));

                centerRV(injectionName) = round(median(runningRV(injectionName).point,1));
                clear dist
                for i=1:size(runningRV(injectionName).point,1)
                    dist(i) = norm(centerRV(injectionName) - runningRV(injectionName).point(i,:));
                end
                meandist = mean(dist);
                %stddist = std(dist);
                stddist = 20;
                flaggedRV = runningRV(injectionName).index(dist>(stddist));
                flagged = [flagged horzcat(flaggedLV,flaggedRV)];
            end
            indexOfFoldersToDo = flagged;
            %FoldersToDo = FoldersToDo(flagged);
            
        end
        if(isempty(flagged) || MyVal == -1)
            MyVal = -1;
        else
            disp('We still have some slices that aren''t consistant with the rest.  Please wait while I redo them');
            usePreviousRVLV = 1;
            askedForParFileOverwrite = 1;
            clear data
            data.bla = 0;
            runningRV = containers.Map({''},{data}); remove(runningRV,'');
            runningLV = containers.Map({''},{data}); remove(runningLV,'');
            MyVal = 1;
        end
    end
    disp('These are the guesses to where the heart is.  The RV point doesn''t matter, while the LV does for some stages.');
    disp('If the cropping is incorrect please change the rangex and rangey in the par file');
    disp('Here are the par files for your convienence');
    myfiles = keys(parFileMap);
    for i=1:length(myfiles)
        disp([myfiles{i} ' : ' parFileMap(myfiles{i})]);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;if ran(stg,0.23,'applies shifts');
    myprocessingDir = pwd;
    shiftfilename = ['Output/shiftsMAN.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.txt'];
    if(~exist(shiftfilename))
        disp('No shift file found');
        return;
    end
    shifts = load(shiftfilename);
    %just do mpi2d.par
    
    dicomfiles = dir([infile '*.dcm']);
    if(~isempty(strfind(infile,'.mat')) || isempty(dicomfiles))
        
        if(exist(['../Archive/Manual/series' num2str(seriesNum) '.slice' num2str(sliceNum) '.par']))
            ParFileName = ['../Archive/Manual/series' num2str(seriesNum) '.slice' num2str(sliceNum) '.par'];
            ReadPar;
        end
        mycine = RetrieveMovieFromMat(infile,rangex, rangey,ranget,flipy,flipx,rotNinety,inendi,intype,inpnxy);
        cinemri1 = zeros(size(mycine));
        for t=1:size(mycine,3)
            temp = circshift(mycine(:,:,t),shifts(t,:));
            cinemri1(:,:,t) = temp;
        end
    else
        cd(infile);
        dicomfiles = dir([infile '*.dcm']);
        clear cinemri1
        % BETTER to take the pains to get recons the size desired... 
        % does make me think though, if Taeho recons frames 2:end, do we
        % put on dicom headers 2:end, or 1:end-1??  Just would affect
        % timestamps
%         if( length(shifts)> length(dicomfiles) )   %EVRD 1/5/11, for recons one shorter. BE CAREFUL with TIMESTAMPS NOW!
%             shifts=shifts(2:end,:); % drop first shifts
%         end
       % for j=1:length(shifts)
        counter=1;
     
        if prod(size(framesToSelect))==0  % trying this EVRD 10/30/11
            framesToSelect=ranget;
        end
        for j=framesToSelect  %ranget    % EVRD 4/26/11, to crop to ranget. Assume subset of shifts. 
            if (j <=max(ranget)) && (j >=min(ranget))
                temp = mpi_upsampleImages(dicomread(dicomfiles(j).name));
               temp = circshift(temp,shifts(counter,:));
               cinemri1(:,:,counter) = temp(rangex, rangey);
                counter=counter+1;
            end
        end
    end
    cd(myprocessingDir);
    save(['Output/cinemri1.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.mat'],'cinemri1');
        
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;if ran(stg,0.26,'make shift file''s reference frame to be [0 0] by subtracting that frame from whole vector');
    myprocessingDir = pwd;
    shiftfilename = ['Output/shiftsMAN.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.txt'];
    if(~exist(shiftfilename))
        disp('No shift file found');
        return;
    end
    shifts = load(shiftfilename);
    if(any(shifts(referenceFrame,:)))
        %compensation
        if(referenceFrame < 0)
            disp('Cannot compensate to Reference frame');
            keyboard;
        end
        for t=1:length(shifts)
            shifts(t,:) = shifts(t,:) - shifts(referenceFrame,:);
        end
        fid = fopen(strcat(outpath,'shiftsMAN.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.txt'),'w');
        for t=1:length(shifts)
            fprintf(fid,'%f	%f\n',shifts(t,:));
        end
        fclose(fid);
        
        dicomfiles = dir([infile '*.dcm']);
        if(~isempty(strfind(infile,'.mat')) || isempty(dicomfiles))

            if(exist(['../Archive/Manual/series' num2str(seriesNum) '.slice' num2str(sliceNum) '.par']))
                ParFileName = ['../Archive/Manual/series' num2str(seriesNum) '.slice' num2str(sliceNum) '.par']; %#ok<NASGU>
                ReadPar;
            end
            mycine = RetrieveMovieFromMat(infile,rangex, rangey,ranget,flipy,flipx,rotNinety,inendi,intype,inpnxy);
            cinemri1 = zeros(size(mycine));
            for t=1:size(mycine,3)
                temp = circshift(mycine(:,:,t),shifts(t,:));
                cinemri1(:,:,t) = temp;
            end
        else
            cd(infile);
            dicomfiles = dir([infile '*.dcm']);
            clear cinemri1
            for j=1:length(shifts)
                temp = mpi_upsampleImages(dicomread(dicomfiles(j).name));
                temp = circshift(temp,shifts(j,:));
                cinemri1(:,:,j) = temp(rangex, rangey);
            end
        end

        cd(myprocessingDir);
        save(['Output/cinemri1.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.mat'],'cinemri1');
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;if ran(stg,0.24,'copy shifts and contours');
    parfiles = dir('*.par');
    previousInfileName = '';
    injections = containers.Map({''},{{''}}); remove(injections,'');
    parfilesBySlice = containers.Map({0},{{''}}); remove(parfilesBySlice,0);
    for i=1:length(parfiles)
        %get at the series number and slice number
        ParFileName = parfiles(i).name; %#ok<NASGU>
        ReadPar;
        
        %group them by slice number
        if(~isKey(parfilesBySlice,sliceNum))
            parfilesBySlice(sliceNum) = {};
        end
        parfileslist = parfilesBySlice(sliceNum);
        parfileslist = vertcat(parfileslist,parfiles(i).name); %#ok<AGROW>
        parfilesBySlice(sliceNum) = parfileslist;
        
        %group them into injections
        genericInfilename = infile(1:strfind(infile,'('));
        if(~strcmp(genericInfilename,previousInfileName))
            previousInfileName = genericInfilename;
            injections(genericInfilename) = {};
        end
        injectionParFiles = injections(genericInfilename);
        injectionParFiles{sliceNum} = parfiles(i).name;
        injections(genericInfilename) = injectionParFiles;
    end
    
    %spread the joy of shifts
    injectionkeys = keys(injections);
    for i=1:length(injectionkeys);
        parfiles = injections(injectionkeys{i});
        goodshiftfilename = '';
        for j=1:length(parfiles)
            if(isempty(parfiles{j})), continue; end;
            ParFileName = parfiles{j}; %#ok<NASGU>
            ReadPar
            shiftfilename = ['Output/shiftsMAN.study' num2str(studyNum) '.slice' num2str(sliceNum) '.txt'];
            if(exist(shiftfilename))
                goodshiftfilename = shiftfilename;
                break;
            end
        end
        if(~strcmp(goodshiftfilename,''))
            %we have a good shift file thus we need top copy it to the
            %other series in this injection
            for j=1:length(parfiles)
                if(isempty(parfiles{j})), continue; end;
                ParFileName = parfiles{j}; %#ok<NASGU>
                ReadPar
                shiftfilename = ['Output/shiftsMAN.study' num2str(studyNum) '.slice' num2str(sliceNum) '.txt'];
                if(~exist(shiftfilename))
                    copyfile(goodshiftfilename,shiftfilename);
                    %now apply those shifts
                    %mpi2d stg=0.23
                end
            end
        end
    end
    
    %go through each slice and look for contours
    slices = keys(parfilesBySlice);
    for i = 1:length(slices)
        slice = slices{i};
        parfileNamesWithSimilarSlice = parfilesBySlice(slice);
        
        copyFromSeriesNum = -1;
        for j=1:length(parfileNamesWithSimilarSlice)
            ParFileName = parfileNamesWithSimilarSlice{j}; %#ok<NASGU>
            ReadPar;
            bloodfile2 = ['blood_polyCoords.study' int2str(studyNum) '.slice' int2str(slice)];
            if(exist(['Output/' bloodfile2]))
                copyFromSeriesNum = studyNum;
                break;
            end
        end
        if(copyFromSeriesNum > 0)
            
            bloodfile = ['blood_polyCoords.study' int2str(copyFromSeriesNum) '.slice' int2str(slice)];
            endofile = ['endo_polyCoords.study' int2str(copyFromSeriesNum) '.slice' int2str(slice)];
            epifile = ['epi_polyCoords.study' int2str(copyFromSeriesNum) '.slice' int2str(slice)];
            Roifile = ['Roi_start_angle.study' int2str(copyFromSeriesNum) '.slice' int2str(slice)];
            for j=1:length(parfileNamesWithSimilarSlice)
                ParFileName = parfileNamesWithSimilarSlice{j}; %#ok<NASGU>
                ReadPar;
                if(studyNum == copyFromSeriesNum), continue; end;
                cd('Output');
                bloodfile2 = ['blood_polyCoords.study' int2str(studyNum) '.slice' int2str(slice)];
                endofile2 = ['endo_polyCoords.study' int2str(studyNum) '.slice' int2str(slice)];
                epifile2 = ['epi_polyCoords.study' int2str(studyNum) '.slice' int2str(slice)];
                Roifile2 = ['Roi_start_angle.study' int2str(studyNum) '.slice' int2str(slice)];
                
                if(~exist(bloodfile2))
                    copyfile(bloodfile,bloodfile2);
                    copyfile(endofile,endofile2);
                    copyfile(epifile,epifile2);
                    copyfile(Roifile,Roifile2);
                end
                cd('..');
            end
        end
    end
    
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;if ran(stg,0.3,'subset pre-Processing');
    SubsetDirectory = '../ReconData/Subsets/';
    %SubsetDirectory = '../ReconData/Subsets/TCR/';
    if(~exist(SubsetDirectory))
        disp('I cannot find the subset directory.  It should be in ReconData/Subsets');
        return;
    end
    FoldersToCheck = dir([SubsetDirectory '/*']);
    ToProcess = containers.Map({''},{{''}});
    ToProcess.remove('');
    mypwd = pwd;
    for i=3:length(FoldersToCheck)
        if(~FoldersToCheck(i).isdir), continue; end;
        if(~strcmp(DoThisFolder,'') && ~strcmp(DoThisFolder,FoldersToCheck(i).name))
            continue;
        end
        dicomFiles = dir([SubsetDirectory '/' FoldersToCheck(i).name '/*.dcm']);
        if(isempty(dicomFiles)) 
            disp([FoldersToCheck(i).name ' contains no Dicom Files.  Skipping']);
            continue; 
        end;
        %add this folder and it's associated dicom files to the map
        ToProcess(FoldersToCheck(i).name) = {dicomFiles(:).name};
    end
    remove(ToProcess,{''});
    
    foldersToDo = keys(ToProcess);
    if(isempty(foldersToDo))
        disp(['I can''t find ' DoThisFolder]);
    end
    for folderi = 1:length(foldersToDo)
        folder = foldersToDo{folderi};
        disp('---------------');
        disp(folder);
        
        %extract the series number and look up mother and slice and copy
        %appropiate files over
        pareni = strfind(folder,'(');
        endpareni = strfind(folder,')');
        series = str2num(folder((pareni+1):(endpareni-1)));
        motherSeries = floor(series/1000)*1000;
        subset = series-motherSeries;
        MotherFolder = strrep(folder,num2str(series),num2str(motherSeries));
        
        if(~exist(['../ReconData/' MotherFolder]))
            disp(['Couldn''t find ' MotherFolder ' in ReconData.  Please put it there so I can figure out what slice this subset it']);
            continue;
        end
        %construct string that will find all slices (share the folder name)
        Regex = strrep(folder,num2str(series),'*');
        candidates = dir(['../ReconData/' Regex]);
        slice = -1;
        for i=1:length(candidates)
            if(strcmp(candidates(i).name,MotherFolder))
                slice = i;
                break;
            end
        end
        if(slice == -1)
            disp(['I could not figure out what slice ' MotherFolder ' is']);
        else
            disp(['I think I am slice# ' num2str(slice)]);
        end
        
        %now on to copying the files
        source = ['.study' num2str(motherSeries) '.slice' num2str(slice)];
        destination = ['.study' num2str(series) '.slice' num2str(slice)];
        
        if(exist(['timeStampSer' num2str(series) '.mat'])), delete(['timeStampSer' num2str(series) '.mat']); end; %#ok<*EXIST>
        copyfile(['timeStampSer' num2str(motherSeries) '.mat'],['timeStampSer' num2str(series) '.mat']);
        cd([mypwd '/Output']);
            shiftfile = ['shiftsMAN' source '.txt'];
            bloodfile = ['blood_polyCoords' source];
            endofile = ['endo_polyCoords' source];
            epifile = ['epi_polyCoords' source];
            roifile = ['Roi_start_angle' source];
            shiftfile2 = ['shiftsMAN' destination '.txt'];
            bloodfile2 = ['blood_polyCoords' destination];
            endofile2 = ['endo_polyCoords' destination];
            epifile2 = ['epi_polyCoords' destination];
            roifile2 = ['Roi_start_angle' destination];
            if(exist(shiftfile2)), delete(shiftfile2); end;
            if(exist(bloodfile2)), delete(bloodfile2); end;
            if(exist(endofile2)), delete(endofile2); end;
            if(exist(epifile2)), delete(epifile2); end;
            if(exist(roifile2)), delete(roifile2); end;
            if(exist(bloodfile))
                copyfile(shiftfile,shiftfile2);
                copyfile(bloodfile,bloodfile2);
                copyfile(endofile,endofile2);
                copyfile(epifile,epifile2);
                copyfile(roifile,roifile2);
            else
                disp(['I couldn''t find the contour files for series:' num2str(series) 'please put it in Processing/Output']);
            end
        cd('..');
       
        
        %now on to creating the cine file
        clear cinemri1 time instance
        dicomFilesToDo = ToProcess(folder);
        %cd(['../ReconData/Subsets/' folder]);
        cd([SubsetDirectory folder]);
        
        disp('Reading in dcm files');
        myinfolder = pwd;
        for dicomi=1:length(dicomFilesToDo)
            cinemri1(:,:,dicomi) = dicomread(dicomFilesToDo{dicomi});
            Header = dicominfo(dicomFilesToDo{dicomi});
            time(dicomi) = str2num(Header.AcquisitionTime);
            instance(dicomi) = Header.InstanceNumber;
        end
        cd(mypwd);
        [instance IX] = sort(instance);
        cinemri1 = cinemri1(:,:,IX);
        time = sort(time);
        
        %now we read in the par file and get the rangex, rangey and do
        %upsampling
        ParFileName = ['series' num2str(motherSeries) '.slice' num2str(slice) '.par']; %#ok<NASGU>
        ReadPar
        
        
        disp(['Upsampling ' folder]);
        cinemri1 = mpi_upsampleImages(cinemri1);
        %grab the shift file and perform the shifts
    %    Shifts = load(['Output/' shiftfile]);  %EVRD removed 12/30/10, not
    %    being used
        %for t=1:length(Shifts)
        %    cinemri1(:,:,t) = circshift(cinemri1(:,:,t),Shifts(t,:));
        %end
        %create the subset par file
        parfile = makeParFile(myinfolder, mypwd, 'Template.par', slice, subset, series,max(ranget),rangex,rangey); %#ok<NASGU>

        %now cut out the heart
        %cinemri1 = cinemri1(rangex,rangey,1:length(Shifts));
        %if(exist(['Output/cinemri1' destination '.mat'])), delete(['Output/cinemri1' destination '.mat']); end;
        %save(['Output/cinemri1' destination '.mat'],'cinemri1');
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end;if ran(stg,1,'INPUT AND SHOW MRI DATA');
	disp('Comment: All future stages depend on this');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end;if ran(stg,1.1,'Show movie of MRI data to select a Reference Frame and to determine cropping dimensions');
% 
%   infile
%   if (infile(length(infile)-3:end) ~= '.mat') %% assumes the file ends with a number, implying it's a DICOM image
%       ff=fopen(infile,'r',inendi);
%       a=fread(ff,intype);
%       fclose(ff);
%       n3=length(a)/prod(inpnxy);
%       a=double(reshape(a,[inpnxy n3]));
%       if(flipx)
%           cinemri=a(max(rangex):-1:min(rangex),rangey,ranget); % this drops out unwanted frames but need to worry about delta_t!!!!
%       else
%           cinemri=a(rangex,rangey,ranget);
%       end
%       if(rotNinety)
%           for ii=1:max(ranget)-min(ranget)+1
%               cinemri(:,:,ii)=rot90(cinemri(:,:,ii),-1);
%           end
%       end
%       
%       for i=1:(max(ranget)-min(ranget)+1)         % transpose for display purposes
%           tmpp(:,:,i)=cinemri(:,:,i)';
%       end
%       cinemri=tmpp;
%       clf;
%       imagesc(cinemri(:,:,referenceFrame)),brighten(.25),axis image,colormap gray,colorbar;
%   else
%       filename=infile; %%% assuming the images are in .mat format and can be loaded directly
%       load (filename)
%       a=imgrr(rangey,rangex,ranget);
%       
%       if(flipy)
%           cinemriy=a(end:-1:1,:,:); % this drops out unwanted frames but need to worry about delta_t!!!!
%       else
%           cinemriy=a(:,:,:);
%       end
%       %%% often the reconstructed .mat images require 'flipping' in the x- and y-directions
%       if(flipx)
%           cinemri3=cinemriy(:,end:-1:1,:); % this drops out unwanted frames but need to worry about delta_t!!!!
%       else
%           cinemri3=cinemriy;
%       end
%       
%       
%       if(rotNinety)
%           for ii=1:max(ranget)-min(ranget)+1
%               cinemri(:,:,ii)=rot90(cinemri3(:,:,ii),-1);
%           end
%       else
%           cinemri(:,:,:)=cinemri3(:,:,:);
%       end
%       
%       for i=1:(max(ranget)-min(ranget)+1)         % transpose for display purposes
%           tmpp(:,:,i)=cinemri(:,:,i);
%       end
%       cinemri=tmpp;
%       clf;
%       imagesc(cinemri(:,:,referenceFrame)),brighten(.25),axis image,colormap gray,colorbar;
%   end
%   
%   title('Reference Frame')
%   disp('Hit key to go on to movie display of frames')
%   pause
%   disp('Choose another reference frame to replace in .par file if desired')
%   showmovie(cinemri,delay); 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
end;if ran(stg,2,'REGISTRATION OF MRI DATA');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end; if ran(stg,2.1,'Perform complete Manual Registration');
    
infile %#ok<NOPRT>
if(exist(infile) == 7)
    infile = [infile '/'];
end
if (infile(end) == '/') %% assumes the file ends with a /, implying it's a folder of DICOM images
    files = dir([infile '*.dcm']);
    [s1 s2] = size(mpi_upsampleImages(dicomread([infile char(files(1).name)])));
    a = zeros(s1,s2,length(files));
    for i=1:length(files)
        a(:,:,i) = mpi_upsampleImages(dicomread([infile char(files(i).name)]));
    end
    if(flipx)
        cinemri=a(max(rangex):-1:min(rangex),rangey,ranget); % this drops out unwanted frames but need to worry about delta_t!!!!
    else
        cinemri=a(rangex,rangey,ranget);
    end
    if(rotNinety)
        for ii=1:max(ranget)-min(ranget)+1
            cinemri(:,:,ii)=rot90(cinemri(:,:,ii),-1);
        end
    end
    
%     for i=1:(max(ranget)-min(ranget)+1)         % transpose for display purposes
%         tmpp(:,:,i)=cinemri(:,:,i)';
%     end
%     cinemri=tmpp;
elseif (infile(length(infile)-3:end) ~= '.mat') %#ok<STCMP> %% assumes the file ends with a number, implying it's a DICOM image
    ff=fopen(infile,'r',inendi);
    a=fread(ff,intype);
    fclose(ff);
    
    n3=length(a)/prod(inpnxy);
    a=double(reshape(a,[inpnxy n3]));
    if(flipx)
        cinemri=a(max(rangex):-1:min(rangex),rangey,ranget); % this drops out unwanted frames but need to worry about delta_t!!!!
    else
        cinemri=a(rangex,rangey,ranget);
    end
    if(rotNinety)
        for ii=1:max(ranget)-min(ranget)+1
            cinemri(:,:,ii)=rot90(cinemri(:,:,ii),-1);
        end
    end
    
    for i=1:(max(ranget)-min(ranget)+1)         % transpose for display purposes
        tmpp(:,:,i)=cinemri(:,:,i)';
    end
    cinemri=tmpp;
else
    filename=infile; %%% assuming the images are in .mat format and can be loaded directly
    load (filename)
    a=imgrr(rangey,rangex,ranget);
    
    if(flipy)
        cinemriy=a(end:-1:1,:,:); % this drops out unwanted frames but need to worry about delta_t!!!!
    else
        cinemriy=a(:,:,:);
    end
    %%% often the reconstructed .mat images require 'flipping' in the x- and y-directions
    if(flipx)
        cinemri3=cinemriy(:,end:-1:1,:); % this drops out unwanted frames but need to worry about delta_t!!!!
    else
        cinemri3=cinemriy;
    end
    
    
    if(rotNinety)
        for ii=1:max(ranget)-min(ranget)+1
            cinemri(:,:,ii)=rot90(cinemri3(:,:,ii),-1);
        end
    else
        cinemri(:,:,:)=cinemri3(:,:,:);
    end
    
    for i=1:(max(ranget)-min(ranget)+1)         % transpose for display purposes
        tmpp(:,:,i)=cinemri(:,:,i);
    end
    cinemri=tmpp;
end

  % added 6/29/04 evrd (This is to upsample to then perform manual registration)
if(flagUpsampleImage==1) 
%  cinemri=mpi_upsampleImages(cinemri);
end

%  cine_up=Supersample(cinemri,up_factor,studyNum,sliceNum);
%  manualShiftsVector=manualRegistration(cinemri,sliceNum, seriesNum,outpath, referenceFrame,1);
manualShiftsVector = auto_registration_trial1(cinemri,ranget,1);

%create a MANshifts file
MANshiftfilename=strcat('shiftsMAN.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.txt');
fid = fopen([outpath '/' MANshiftfilename],'w');
for t=1:length(MANshiftfilename)
    fprintf(fid,'%f %f\n',manualShiftsVector(t,1),manualShiftsVector(t,2));
end
fclose(fid);
 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end;if ran(stg,2.5,'Mutual Information');

    disp('Retrieving images and upsampling...');
    dicoms = dir([infile '*.dcm']);
    legacy = 0;
    if(isempty(dicoms))
        %Now check if this is a legacy par file
        if(exist(infile))
            legacy = 1;
            %this is probably either a mat or a raw file
            cinemri1 = RetrieveMovieFromMat(infile,rangex, rangey,ranget,flipy,flipx,rotNinety,inendi,intype,inpnxy);
            mycine = cinemri1;
            [sx sy ~] = size(cinemri1);
            if(exist(['../Archive/Manual/series' num2str(seriesNum) '.slice' num2str(sliceNum) '.par']))
                ParFileName = ['../Archive/Manual/series' num2str(seriesNum) '.slice' num2str(sliceNum) '.par'];
                ReadPar;
            end
            oldrangex = rangex;
            oldrangey = rangey;
            ParFileName = 'mpi2d.par';
            ReadPar
            rangex = 1:sx;
            rangey = 1:sy;
        else
            disp('No dicoms found.  Cannot do cross correlation registration');
            keyboard;
            return;
        end
    else
        for t=1:length(dicoms)
            mycine(:,:,t) = mpi_upsampleImages(dicomread([infile dicoms(t).name]));
            HeaderInfo = dicominfo([infile dicoms(t).name]);
            mytime(t) = str2num(HeaderInfo.AcquisitionTime); %#ok<ST2NM>
        end
        [mytime, IX] = sort(mytime);
        mycine = mycine(:,:,IX);
        if(referenceFrame < 0)
            if(exist(['Output/RVLVpos.study' int2str(seriesNum) '.slice' int2str(sliceNum) '.mat']))
                load(['Output/RVLVpos.study' int2str(seriesNum) '.slice' int2str(sliceNum) '.mat']);
            else
                [RV,LV] = FindLVRV(mycine,0);
            end
            if(referenceFrame == -1)
                [~, referenceFrame] = max(mycine(RV(1), RV(2),:).*mycine(LV(1), LV(2),:));
            elseif(referenceFrame == -2)
                [~, referenceFrame] = max(mycine(LV(1), LV(2),:));
            else
                referenceFrame = size(mycine,3);
            end
        end
    end
    
	shifts_out = mutualInformation(mycine, rangex, rangey,referenceFrame);
    
    if(legacy)
        %nate doesn't choose a rangex that's narrow enough for most
        %registration algorithms.  We must change it so that the algorithm
        %does well.  This is here so we create a cine file that makes the
        %contours overlay the heart properly.
        rangex = oldrangex;
        rangey = oldrangey;
        cinemri1 = RetrieveMovieFromMat(infile,rangex, rangey,ranget,flipy,flipx,rotNinety,inendi,intype,inpnxy);
        mycine = cinemri1;
    end
    clear cinemri1
    for t=1:length(shifts_out)
        temp = circshift(mycine(:,:,t),shifts_out(t,:));
        if(legacy)
            cinemri1(:,:,t) = temp(:, :);
        else
            cinemri1(:,:,t) = temp(rangex, rangey);
        end
    end
    disp('Saving cinemri file');
    save(['Output/cinemri1.study' int2str(seriesNum) '.slice' int2str(sliceNum) '.mat'],'cinemri1');
    disp('Saving shift file');
    fid = fopen(strcat(outpath,'shiftsMAN.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.txt'),'w');
    for t=1:length(shifts_out)
        fprintf(fid,'%f	%f\n',shifts_out(t,:));
    end
    fclose(fid);
    figure(seriesNum*10 + 2*sliceNum), subplot(2,2,2), cla, hold on;
    plot(shifts_out(:,1));
    plot(shifts_out(:,2)+10,'g');
    title('Shifts(Bottom:X, Top:Y)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end;if ran(stg,2.6,'xcoor(herusitic)');

    
    disp('Retrieving images and upsampling...');
    dicoms = dir([infile '*.dcm']);
    legacy = 0;
    if(isempty(dicoms))
        %Now check if this is a legacy par file
        if(exist(infile))
            legacy = 1;
            %this is probably either a mat or a raw file
            cinemri1 = RetrieveMovieFromMat(infile,rangex, rangey,ranget,flipy,flipx,rotNinety,inendi,intype,inpnxy);
            mycine = cinemri1;
            [sx sy ~] = size(cinemri1);
            if(exist(['../Archive/Manual/series' num2str(seriesNum) '.slice' num2str(sliceNum) '.par']))
                ParFileName = ['../Archive/Manual/series' num2str(seriesNum) '.slice' num2str(sliceNum) '.par'];
                ReadPar;
            end
            oldrangex = rangex;
            oldrangey = rangey;
            ParFileName = 'mpi2d.par';
            ReadPar
            rangex = 1:sx;
            rangey = 1:sy;
        else
            disp('No dicoms found.  Cannot do cross correlation registration');
            keyboard;
            return;
        end
    else
        for t=1:length(dicoms)
            mycine(:,:,t) = mpi_upsampleImages(dicomread([infile dicoms(t).name]));
            HeaderInfo = dicominfo([infile dicoms(t).name]);
            mytime(t) = str2num(HeaderInfo.AcquisitionTime); %#ok<ST2NM>
        end
        [mytime, IX] = sort(mytime);
        mycine = mycine(:,:,IX);
        if(referenceFrame < 0)
            if(exist(['Output/RVLVpos.study' int2str(seriesNum) '.slice' int2str(sliceNum) '.mat']))
                load(['Output/RVLVpos.study' int2str(seriesNum) '.slice' int2str(sliceNum) '.mat']);
            else
                [RV,LV] = FindLVRV(mycine,0);
            end
            if(referenceFrame == -1)
                [~, referenceFrame] = max(mycine(RV(1), RV(2),:).*mycine(LV(1), LV(2),:));
            elseif(referenceFrame == -2)
                [~, referenceFrame] = max(mycine(LV(1), LV(2),:));
            else
                referenceFrame = size(mycine,3);
            end
        end
    end
    
    
    %do xcoor here
    corr_offset(size(mycine,3),:) = [0 0];
    [sx sy st] = size(mycine(rangex, rangey,:));
    %create filter
    reference = mycine(rangex,rangey,referenceFrame);
    c = normxcorr2(mycine(rangex,rangey,4),reference);
    
    for t=(size(mycine,3)-1):-1:1
        %referenceFrame = mycine(rangex,rangey,t+1);
        reference = mycine(rangex,rangey,referenceFrame);
        c = normxcorr2(mycine(rangex,rangey,t),reference);
        h = zeros(size(c));
        width = round(min(t/2,30));
        range1 = (sx-width):(sx+width);
        range2 = (sy-width):(sy+width);
        h(range1,range2) = 1;
        c = c .* h;
        [~, imax] = max(abs(c(:)));
        [ypeak, xpeak] = ind2sub(size(c),imax(1));
        corr_offset(t,:) = [(ypeak-sy) (xpeak-sx)];
        if(norm(corr_offset(t,:)) > 20)
            corr_offset(t,:) = [0 0];
        end
        mycine(:,:,t) = circshift(mycine(:,:,t), corr_offset(t,:));
    end
    shifts_out = corr_offset;
    manhatten = sum(abs(shifts_out - repmat(median(shifts_out,1),size(mycine,3),1)),2);
    if(length(manhatten(manhatten>0)) > length(manhatten)/3)
        for t=1:st
            if(manhatten(t)>15)
                shifts_out(t,:) = [0 0];
            end
        end
    end
    
    
    
    if(legacy)
        %nate doesn't choose a rangex that's narrow enough for most
        %registration algorithms.  We must change it so that the algorithm
        %does well.  This is here so we create a cine file that makes the
        %contours overlay the heart properly.
        rangex = oldrangex;
        rangey = oldrangey;
        cinemri1 = RetrieveMovieFromMat(infile,rangex, rangey,ranget,flipy,flipx,rotNinety,inendi,intype,inpnxy);
        mycine = cinemri1;
    end
    clear cinemri1
    for t=1:length(shifts_out)
        temp = circshift(mycine(:,:,t),shifts_out(t,:));
        if(legacy)
            cinemri1(:,:,t) = temp(:, :);
        else
            cinemri1(:,:,t) = temp(rangex, rangey);
        end
    end
    disp('Saving cinemri file');
    save(['Output/cinemri1.study' int2str(seriesNum) '.slice' int2str(sliceNum) '.mat'],'cinemri1');
    
    disp('Saving shift file');
    fid = fopen(strcat(outpath,'shiftsMAN.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.txt'),'w');
    for t=1:length(shifts_out)
        fprintf(fid,'%f	%f\n',shifts_out(t,:));
    end
    fclose(fid);

    figure(seriesNum*10 + 2*sliceNum), subplot(2,2,2), cla, hold on;
    plot(shifts_out(:,1));
    plot(shifts_out(:,2)+10,'g');
    title('Shifts(Bottom:X, Top:Y)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end;if ran(stg,2.61,'xcoor');

    
    disp('Retrieving images and upsampling...');
    dicoms = dir([infile '*.dcm']);
    legacy = 0;
    if(isempty(dicoms))
        %Now check if this is a legacy par file
        if(exist(infile))
            legacy = 1;
            %this is probably either a mat or a raw file
            cinemri1 = RetrieveMovieFromMat(infile,rangex, rangey,ranget,flipy,flipx,rotNinety,inendi,intype,inpnxy);
            mycine = cinemri1;
            [sx sy ~] = size(cinemri1);
            if(exist(['../Archive/Manual/series' num2str(seriesNum) '.slice' num2str(sliceNum) '.par']))
                ParFileName = ['../Archive/Manual/series' num2str(seriesNum) '.slice' num2str(sliceNum) '.par'];
                ReadPar;
            end
            oldrangex = rangex;
            oldrangey = rangey;
            ParFileName = 'mpi2d.par';
            ReadPar
            rangex = 1:sx;
            rangey = 1:sy;
        else
            disp('No dicoms found.  Cannot do cross correlation registration');
            keyboard;
            return;
        end
    else
        for t=1:length(dicoms)
            mycine(:,:,t) = mpi_upsampleImages(dicomread([infile dicoms(t).name]));
            HeaderInfo = dicominfo([infile dicoms(t).name]);
            mytime(t) = str2num(HeaderInfo.AcquisitionTime); %#ok<ST2NM>
        end
        [mytime, IX] = sort(mytime);
        mycine = mycine(:,:,IX);
        if(referenceFrame < 0)
            if(exist(['Output/RVLVpos.study' int2str(seriesNum) '.slice' int2str(sliceNum) '.mat']))
                load(['Output/RVLVpos.study' int2str(seriesNum) '.slice' int2str(sliceNum) '.mat']);
            else
                [RV,LV] = FindLVRV(mycine,0);
            end
            if(referenceFrame == -1)
                [~, referenceFrame] = max(mycine(RV(1), RV(2),:).*mycine(LV(1), LV(2),:));
            elseif(referenceFrame == -2)
                [~, referenceFrame] = max(mycine(LV(1), LV(2),:));
            else
                referenceFrame = size(mycine,3);
            end
        end
    end
    
    
    %do xcoor here
    corr_offset(size(mycine,3),:) = [0 0];
    [sx sy st] = size(mycine(rangex, rangey,:));
    %create filter
    reference = mycine(rangex,rangey,referenceFrame);
    c = normxcorr2(mycine(rangex,rangey,4),reference);
    for t=(size(mycine,3)-1):-1:1
        %referenceFrame = mycine(rangex,rangey,t+1);
        reference = mycine(rangex,rangey,referenceFrame);
        c = normxcorr2(mycine(rangex,rangey,t),reference);
        [~, imax] = max(abs(c(:)));
        [ypeak, xpeak] = ind2sub(size(c),imax(1));
        corr_offset(t,:) = [(ypeak-sy) (xpeak-sx)];
        mycine(:,:,t) = circshift(mycine(:,:,t), corr_offset(t,:));
    end
    shifts_out = corr_offset;
    
    
    if(legacy)
        %nate doesn't choose a rangex that's narrow enough for most
        %registration algorithms.  We must change it so that the algorithm
        %does well.  This is here so we create a cine file that makes the
        %contours overlay the heart properly.
        rangex = oldrangex;
        rangey = oldrangey;
        cinemri1 = RetrieveMovieFromMat(infile,rangex, rangey,ranget,flipy,flipx,rotNinety,inendi,intype,inpnxy);
        mycine = cinemri1;
    end
    clear cinemri1
    for t=1:length(shifts_out)
        temp = circshift(mycine(:,:,t),shifts_out(t,:));
        if(legacy)
            cinemri1(:,:,t) = temp(:, :);
        else
            cinemri1(:,:,t) = temp(rangex, rangey);
        end
    end
    disp('Saving cinemri file');
    save(['Output/cinemri1.study' int2str(seriesNum) '.slice' int2str(sliceNum) '.mat'],'cinemri1');
    
    disp('Saving shift file');
    fid = fopen(strcat(outpath,'shiftsMAN.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.txt'),'w');
    for t=1:length(shifts_out)
        fprintf(fid,'%f	%f\n',shifts_out(t,:));
    end
    fclose(fid);

    figure(seriesNum*10 + 2*sliceNum), subplot(2,2,2), cla, hold on;
    plot(shifts_out(:,1));
    plot(shifts_out(:,2)+10,'g');
    title('Shifts(Bottom:X, Top:Y)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end;if ran(stg,2.8,'Temporal Smoothing');
disp('Retrieving images and upsampling...');
    dicoms = dir([infile '*.dcm']);
    legacy = 0;
    if(isempty(dicoms))
        %Now check if this is a legacy par file
        if(exist(infile))
            legacy = 1;
            %this is probably either a mat or a raw file
            cinemri1 = RetrieveMovieFromMat(infile,rangex, rangey,ranget,flipy,flipx,rotNinety,inendi,intype,inpnxy);
            mycine = cinemri1;
            [sx sy ~] = size(cinemri1);
            if(exist(['../Archive/Manual/series' num2str(seriesNum) '.slice' num2str(sliceNum) '.par']))
                ParFileName = ['../Archive/Manual/series' num2str(seriesNum) '.slice' num2str(sliceNum) '.par'];
                ReadPar;
            end
            oldrangex = rangex;
            oldrangey = rangey;
            ParFileName = 'mpi2d.par';
            ReadPar
            rangex = 1:sx;
            rangey = 1:sy;
        else
            disp('No dicoms found.  Cannot do cross correlation registration');
            keyboard;
            return;
        end
    else
        counter=1;
        if ~isempty(framesToSelect)
          for j=framesToSelect    % EVRD 8/17/11, to select ranget. 
            if (j <=max(ranget)) && (j >=min(ranget))
                temp = mpi_upsampleImages(dicomread([infile dicoms(j).name]));
                HeaderInfo = dicominfo([infile dicoms(j).name]);
                mytime(counter) = str2num(HeaderInfo.AcquisitionTime);
                mycine(:,:,counter) = temp;  % see note just below(rangex, rangey);
                counter=counter+1;
            end
          end
        else
             for j=ranget
                temp = mpi_upsampleImages(dicomread([infile dicoms(j).name]));
                HeaderInfo = dicominfo([infile dicoms(j).name]);
                mytime(counter) = str2num(HeaderInfo.AcquisitionTime);
                mycine(:,:,counter) = temp;  % note need outside crops for Neighboring... (rangex, rangey);
                counter=counter+1;
             end
        end
        [mytime, IX] = sort(mytime); % should we flag if any out of order??
        mycine = mycine(:,:,IX);
        
 
%         for t=1:length(dicoms)
%             mycine(:,:,t) = mpi_upsampleImages(dicomread([infile dicoms(t).name]));
%             HeaderInfo = dicominfo([infile dicoms(t).name]);
%             mytime(t) = str2num(HeaderInfo.AcquisitionTime); %#ok<ST2NM>
%         end
%         [mytime, IX] = sort(mytime);
%         mycine = mycine(:,:,IX);
        if(referenceFrame < 0)
            if(exist(['Output/RVLVpos.study' int2str(seriesNum) '.slice' int2str(sliceNum) '.mat']))
                load(['Output/RVLVpos.study' int2str(seriesNum) '.slice' int2str(sliceNum) '.mat']);
            else
                [RV,LV] = FindLVRV(mycine,0);
            end
            if(referenceFrame == -1)
                [~, referenceFrame] = max(mycine(RV(1), RV(2),:).*mycine(LV(1), LV(2),:));
            elseif(referenceFrame == -2)
                [~, referenceFrame] = max(mycine(LV(1), LV(2),:));
            else
                referenceFrame = size(mycine,3);
            end
        end
    end
    
    %x = load(['timeStampSer' num2str(seriesNum) '.mat']);
    x = load(timeStampFile);  %EVRD 8/11
    timestamps = x.timeStamp;
    if(length(timestamps) ~= size(mycine,3))
        timestamps = timestamps(min(ranget):end);
        if(length(timestamps) ~= size(mycine,3))
            timestamps = x.timeStamp(ranget);
            if(length(timestamps) ~= size(mycine,3))
                disp('Houston, we have a timestamp misalignment problem - can skip, doesnt use in 2.8');
                %keyboard;
            end
        end
    end
    
	shifts_out = NeighboringTemporalSmoothingRegistration(mycine, rangex, rangey, 10, timestamps); % doesn't use timestamps
 
    if(legacy)
        %nate doesn't choose a rangex that's narrow enough for most
        %registration algorithms.  We must change it so that the algorithm
        %does well.  This is here so we create a cine file that makes the
        %contours overlay the heart properly.
        rangex = oldrangex;
        rangey = oldrangey;
        cinemri1 = RetrieveMovieFromMat(infile,rangex, rangey,ranget,flipy,flipx,rotNinety,inendi,intype,inpnxy);
        mycine = cinemri1;
    end
    clear cinemri1
    for t=1:length(shifts_out)
        temp = circshift(mycine(:,:,t),shifts_out(t,:));
        if(legacy)
            cinemri1(:,:,t) = temp(:, :);
        else
            cinemri1(:,:,t) = temp(rangex, rangey);
        end
    end
    disp('Saving cinemri file');
    save(['Output/cinemri1.study' int2str(seriesNum) '.slice' int2str(sliceNum) '.mat'],'cinemri1');
    disp('Saving shift file');
    fid = fopen(strcat(outpath,'shiftsMAN.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.txt'),'w');
    for t=1:length(shifts_out)
        fprintf(fid,'%f	%f\n',shifts_out(t,:));
    end
    fclose(fid);

    figure(seriesNum*10 + 2*sliceNum), subplot(2,2,2), cla, hold on;
    plot(shifts_out(:,1));
    plot(shifts_out(:,2)+10,'g');
    title('Shifts(Bottom:X, Top:Y)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end;if ran(stg,2.9,'Model-based registration');

    
    disp('Retrieving images and upsampling...');
    if(~isempty(strfind(infile,'.mat')))
        mycine = RetrieveMovieFromMat(infile,rangex, rangey,ranget,flipy,flipx,rotNinety,inendi,intype,inpnxy);
        for t=1:size(mycine,3)
            fullcinemri1(:,:,t) = mpi_upsampleImages(mycine(:,:,t));
        end
    else
        dicoms = dir([infile '*.dcm']);
        if(isempty(dicoms))
            disp('No dicoms found.  Cannot do Temporal Smoothing registration');
            keyboard;
            return;
        end
        for t=1:length(dicoms)
            fullcinemri1(:,:,t) = mpi_upsampleImages(dicomread([infile dicoms(t).name]));
            HeaderInfo = dicominfo([infile dicoms(t).name]);
            mytime(t) = str2num(HeaderInfo.AcquisitionTime);
        end
        [mytime, IX] = sort(mytime); %#ok<ASGLU>
        fullcinemri1 = fullcinemri1(:,:,IX);
    end
    
    %% modelBased Registration
    figure(seriesNum*10 + 2*sliceNum)
    %[shifts,parameters,delays,pixelwisevariance,cinemri1] = modelRegistration(fullcinemri1,rangex, rangey,0,2); %#ok<ASGLU>
    [shifts,pixelwisevariance] = model_newer(fullcinemri1, rangex, rangey, 8);
    [~,~,st] = size(cinemri1);

    figure(seriesNum*10 + 2*sliceNum),h = subplot(2,2,3); imagesc(delays), colormap hsv;colorbar; title('pixelwise Delays(time @ onset)');cbfreeze(h);
    mycine = zeros(length(rangex),length(rangey),st);
    for t=1:length(shifts)
        temp = circshift(fullcinemri1(:,:,t),shifts(t,:));
        mycine(:,:,t) = temp(rangex,rangey);
    end
    cinemri1 = mycine;
    
    save(['Output/cinemri1.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.mat'],'cinemri1');
    fid = fopen(['Output/shiftsMAN.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.txt'],'w');
    for t=1:length(shifts)
        fprintf(fid,'%f	%f\n',shifts(t,:));
    end
    fclose(fid);

    figure(seriesNum*10 + 2*sliceNum), subplot(2,2,2), cla, hold on;
    plot(shifts(:,1));
    plot(shifts(:,2)+10,'g');
    title('Shifts(Bottom:X, Top:Y)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end;if ran(stg,2.91,'Ganesh Model-based registration');

    
    disp('Retrieving images and upsampling...');
    if(~isempty(strfind(infile,'.mat')))
        mycine = RetrieveMovieFromMat(infile,rangex, rangey,ranget,flipy,flipx,rotNinety,inendi,intype,inpnxy);
        for t=1:size(mycine,3)
            fullcinemri1(:,:,t) = mpi_upsampleImages(mycine(:,:,t));
        end
    else
        dicoms = dir([infile '*.dcm']);
        if(isempty(dicoms))
            disp('No dicoms found.  Cannot do Temporal Smoothing registration');
            keyboard;
            return;
        end
        for t=1:length(dicoms)
            fullcinemri1(:,:,t) = mpi_upsampleImages(dicomread([infile dicoms(t).name]));
            HeaderInfo = dicominfo([infile dicoms(t).name]);
            mytime(t) = str2num(HeaderInfo.AcquisitionTime);
        end
        [mytime, IX] = sort(mytime); %#ok<ASGLU>
        fullcinemri1 = fullcinemri1(:,:,IX);
        keyboard 
    end
    
    %% modelBased Registration
    figure(seriesNum*10 + 2*sliceNum)
    if(exist(['Output/RVLVpos.study' int2str(seriesNum) '.slice' int2str(sliceNum) '.mat']))
        load(['Output/RVLVpos.study' int2str(seriesNum) '.slice' int2str(sliceNum) '.mat']);
    else
        [RV,LV] = FindLVRV(fullcinemri1(rangex,rangey,:),0);
    end
    mpi2d stg=2.8
    shifts = load(['Output/shiftsMAN.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.txt']);
    
    for t=1:length(shifts)
        fullcinemri2(:,:,t) = circshift(fullcinemri1(:,:,t),shifts(t,:));
    end
    
    shifts2 = GaneshModel(fullcinemri2, rangex, rangey, ranget,RV,referenceFrame);
    
    shifts = shifts(ranget,:) + shifts2;
    [~,~,st] = size(fullcinemri1);

    mycine = zeros(length(rangex),length(rangey),st);
    for t=1:length(shifts)
        temp = circshift(fullcinemri1(:,:,t),shifts(t,:));
        mycine(:,:,t) = temp(rangex,rangey);
    end
    cinemri1 = mycine;
    
    save(['Output/cinemri1.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.mat'],'cinemri1');
    fid = fopen(['Output/shiftsMAN.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.txt'],'w');
    for t=1:length(shifts)
        fprintf(fid,'%f	%f\n',shifts(t,:));
    end
    fclose(fid);

    figure(seriesNum*10 + 2*sliceNum), subplot(2,2,2), cla, hold on;
    plot(shifts(:,1));
    plot(shifts(:,2)+10,'g');
    title('Shifts(Bottom:X, Top:Y+10)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end;if ran(stg,2.11,'zero-Shift registration');

    
    disp('Retrieving images and upsampling...');
    if(~isempty(strfind(infile,'.mat')))
        mycine = RetrieveMovieFromMat(infile,rangex, rangey,ranget,flipy,flipx,rotNinety,inendi,intype,inpnxy);
        disp('not supported right now, ask ED ')
        keyboard
        for t=1:size(mycine,3)
            fullcinemri1(:,:,t) = mpi_upsampleImages(mycine(:,:,t));
        end
    else
        dicoms = dir([infile '*.dcm']);
        if(isempty(dicoms))
            disp('No dicoms found.  Cannot do zero-shift registration');
            keyboard;
            return;
        end
        counter=1;
        if ~isempty(framesToSelect)
          for j=framesToSelect    % EVRD 8/17/11, to select ranget. 
            if (j <=max(ranget)) && (j >=min(ranget))
                temp = mpi_upsampleImages(dicomread([infile dicoms(j).name]));
                HeaderInfo = dicominfo([infile dicoms(j).name]);
                mytime(counter) = str2num(HeaderInfo.AcquisitionTime);
                cinemri1(:,:,counter) = temp(rangex, rangey);
                counter=counter+1;
            end
          end
        else
             for j=ranget
                temp = mpi_upsampleImages(dicomread([infile dicoms(j).name]));
                HeaderInfo = dicominfo([infile dicoms(j).name]);
                mytime(j) = str2num(HeaderInfo.AcquisitionTime);
                cinemri1(:,:,counter) = temp(rangex, rangey);
                counter=counter+1;
             end
        end
        [mytime, IX] = sort(mytime); % should we flag if any out of order??
        cinemri1 = cinemri1(:,:,IX);
            
        
%         
%         for t=1:length(dicoms)
%             fullcinemri1(:,:,t) = mpi_upsampleImages(dicomread([infile dicoms(t).name]));
%             HeaderInfo = dicominfo([infile dicoms(t).name]);
%             mytime(t) = str2num(HeaderInfo.AcquisitionTime);
%         end
%         [mytime, IX] = sort(mytime); %#ok<ASGLU>
%         fullcinemri1 = fullcinemri1(:,:,IX);
    end
    
    shifts = zeros(size(cinemri1,3),2);
    %cinemri1 = fullcinemri1(rangex, rangey,:);

    save(['Output/cinemri1.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.mat'],'cinemri1');
    fid = fopen(['Output/shiftsMAN.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.txt'],'w');
    for t=1:length(shifts)
        fprintf(fid,'%f	%f\n',shifts(t,:));
    end
    fclose(fid);

    figure(seriesNum*10 + 2*sliceNum), subplot(2,2,2), cla, hold on;
    plot(shifts(:,1));
    plot(shifts(:,2)+10,'g');
    title('Shifts(Bottom:X, Top:Y)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end;if ran(stg,3,'SEGMENTATION OF MRI DATA');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;if ran(stg,3.01,'Draw endocardial, epicardial, and LV blood contours for all slices and injections');
    displayAllMoviesAndContours();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;if ran(stg,3.1,'Draw endocardial, epicardial, and LV blood contours'); %#ok<*ALIGN>

  infilename=strcat(outpath,'cinemri1.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat');
    
  load(infilename);

  mpi_drawEndoAndEpi(cinemri1,sliceNum, seriesNum,outpath, referenceFrame)   % interactive contour drawing, change mpi_roipoly for mouse interaction choice

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;if ran(stg,3.11,'SemiAuto-Segmentation endocardial, epicardial, and LV blood contours');

  %infilename=strcat(outpath,'cinemri1.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat');
    
  handle = CR_editor(seriesNum,sliceNum);
%  auto_segmentation(infilename);
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;if ran(stg,3.12,'Auto-Segmentation endocardial, epicardial, and LV blood contours');

    [blood ,epi,endo,angle] = AutoSegment(seriesNum,sliceNum,myothreshold);
    if(isempty(blood))
        disp('blood Contour not created, should still be able to do 3.11');
        return;
    end
    
    fid = fopen([outpath 'blood_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)],'w');
    if(fid < 0)
      disp('Could not open Blood file for saving contour');
    end
    for i=1:length(blood)
       fprintf(fid,'%d %d\n',blood(i,2:-1:1));
    end
    fclose(fid);

    fid = fopen([outpath 'epi_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)],'w');
    if(fid < 0)
      disp('Could not open epi file for saving contour');
    end
    for i=1:length(epi)
       fprintf(fid,'%d %d\n',epi(i,2:-1:1));
    end
    fclose(fid);

    fid = fopen([outpath 'endo_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)],'w');
    if(fid < 0)
      disp('Could not open endo file for saving contour');
    end
    for i=1:length(endo)
       fprintf(fid,'%d %d\n',endo(i,2:-1:1));
    end
    fclose(fid);

    fid = fopen([outpath 'Roi_start_angle.study' int2str(seriesNum) '.slice' int2str(sliceNum)],'w');
    if(fid < 0)
      disp('Could not open Roi_start_angle file for saving contour');
    end
    fprintf(fid,'%d\n',angle);
    fclose(fid);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;if ran(stg,3.13,'Auto-Segmentation endocardial, epicardial, and LV blood contours with reuse');

    [blood ,epi,endo,angle] = AutoSegment(seriesNum,sliceNum,myothreshold,1);
    if(isempty(blood))
        disp('Contours not created');
        return;
    end
    
    fid = fopen([outpath 'blood_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)],'w');
    if(fid < 0)
      disp('Could not open Blood file for saving contour');
    end
    for i=1:length(blood)
       fprintf(fid,'%d %d\n',blood(i,2:-1:1));
    end
    fclose(fid);

    fid = fopen([outpath 'epi_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)],'w');
    if(fid < 0)
      disp('Could not open epi file for saving contour');
    end
    for i=1:length(epi)
       fprintf(fid,'%d %d\n',epi(i,2:-1:1));
    end
    fclose(fid);

    fid = fopen([outpath 'endo_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)],'w');
    if(fid < 0)
      disp('Could not open endo file for saving contour');
    end
    for i=1:length(endo)
       fprintf(fid,'%d %d\n',endo(i,2:-1:1));
    end
    fclose(fid);

    fid = fopen([outpath 'Roi_start_angle.study' int2str(seriesNum) '.slice' int2str(sliceNum)],'w');
    if(fid < 0)
      disp('Could not open Roi_start_angle file for saving contour');
    end
    fprintf(fid,'%d\n',angle);
    fclose(fid);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;if ran(stg,3.2,'Read in segmentation contours to generate regional SI vs. Time curves');
% Read in contours and display, also creates time curves
  infilename=strcat(outpath,'cinemri1.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat');
  load(infilename);

    % moved to here 10/24/05   give input of vector framesToSkip.
    % Must be numbers relative to ranget  (ranget needs to start at 1...)
    if framesToSkip
      framesToSkip 
      cinemri1=mpi_skipFrames(framesToSkip, cinemri1, ranget);
    end
    
    %   curves=mpi_applyRegions(cinemri1,sliceNum,seriesNum,outpath, referenceFrame, numAzimuthalRegions, flagPixelwise, numRadialRegions, pixelSizeRatio);  
    blood = load([outpath 'blood_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)]);
    endo = load([outpath 'endo_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)]);
    epi = load([outpath 'epi_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)]);
    angle = load([outpath 'Roi_start_angle.study' int2str(seriesNum) '.slice' int2str(sliceNum)]);
    figure(seriesNum*10 + 2*sliceNum+1);
    
    curves = ExtractCurves(cinemri1,blood,endo,epi,angle,numAzimuthalRegions,referenceFrame);
    
    curvefilename=strcat(outpath,'curves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat');
    save(curvefilename, 'curves')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;if ran(stg,3.22,'Read in segmentation contours to generate regional SI vs. Time curves and blind estimate of AIF');
% Read in contours and display, also creates time curves
  infilename=strcat(outpath,'cinemri1.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat');
  load(infilename);

% moved to here 10/24/05   give input of vector framesToSkip.
% Must be numbers relative to ranget  (ranget needs to start at 1...)
if framesToSkip
  framesToSkip %#ok<NOPRT>
  cinemri1=mpi_skipFrames(framesToSkip, cinemri1);
end

%     curves=mpi_applyRegions(cinemri1,sliceNum,seriesNum,outpath, referenceFrame, numAzimuthalRegions, flagPixelwise, numRadialRegions, pixelSizeRatio);
    deltaSIcurves=mpi_si2deltaSI(sliceNum,seriesNum,outpath,numPreContrast_bld, numPreContrast_tiss, numSkip, 0);
    
    LV = load([outpath 'blood_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)]);   
    epi = load([outpath 'epi_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)]); 
    endo = load([outpath 'endo_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)]);
    motherSeries = floor(seriesNum/1000)*1000;   %EVRD 10/9/11
    load(['timeStampSer' num2str(motherSeries) '.mat']);
    close all
    aif = CAMM(cinemri1,{epi,endo,LV},timeStamp);
    if(size(deltaSIcurves,2) ~= length(aif))
        disp('blind estimate gave an aif array longer than that of the tissue curves');
        keyboard;
    end
    deltaSIcurves(1,:) = aif;
   
  curvefilename=strcat(outpath,'deltaSIcurves.study',int2str(seriesNum),'.slice',int2str(sliceNum+20),'.mat');
  save(curvefilename, 'deltaSIcurves')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;if ran(stg,3.21,'Perform a frame-by-frame check of Registration and Segmentation');
% Read in contours and display - still working on improvements
  infilename=strcat(outpath,'cinemri1.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat');
  load(infilename);
  mpi_showRegions_old(cinemri1,sliceNum, seriesNum, outpath, numAzimuthalRegions, numRadialRegions, referenceFrame);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;if ran(stg,3.3,'Display the regional SI vs. Time curves');

  curvefilename=strcat(outpath,'curves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat');
  load(curvefilename);
  showcurves(curves,'Sig. Intensity (arb. units)');     % y label

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;if ran(stg,3.4,'Convert the SI curves to "normalized" deltaSI curves and display');
% convert curves from si to gd conc.
  scaleToStudyOne=0;    % need to put with parameters in beginning. This may just set to zero because some studeis not being processed conventionally. Need to worry about this!
  if  useDeltaSI~=1    %  12/30/03.
     gdcurves=mpi_si2gdconc(sliceNum,seriesNum,outpath,numPreContrast_bld, numPreContrast_tiss, numSkip, scaleToStudyOne); %#ok<NASGU>
     % last param is scale to study1 - - will usually want on (=1)!
     curvefilename=strcat(outpath,'gdcurves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat');
     save(curvefilename, 'gdcurves')
  else
    h = figure;
    deltaSIcurves=mpi_si2deltaSI(sliceNum,seriesNum,outpath,numPreContrast_bld, numPreContrast_tiss, numSkip, scaleToStudyOne);
    close(h);
%     deltaSIcurves=deltaSIcurves(:,ranget); % added 10/6/05 so no need to run 1.1. again
    curvefilename=strcat(outpath,'deltaSIcurves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat');
    save(curvefilename, 'deltaSIcurves');
    mycolors = lines(numAzimuthalRegions);
    figure(seriesNum*10 + 2*sliceNum+1), clf, hold on, plot(deltaSIcurves(1,:),'k');
    for region=1:numAzimuthalRegions
       plot(deltaSIcurves(region+1,:)','Color',mycolors(region,:));
    end
  end
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;if ran(stg,3.41,'Convert the SI curves to "normalized" deltaSI curves and display AFTER Proton Density (PD) correction has been performed');
% convert curves from si to gd conc.
  scaleToStudyOne=0;    % need to put with parameters in beginning. This may just set to zero because some studeis not being processed conventionally. Need to worry about this!
  if  useDeltaSI~=1    %  12/30/03.
     gdcurves=mpi_si2gdconc(sliceNum,seriesNum,outpath,numPreContrast_bld, numPreContrast_tiss, numSkip, scaleToStudyOne); %#ok<NASGU>
     % last param is scale to study1 - - will usually want on (=1)!
     curvefilename=strcat(outpath,'gdcurves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat');
     save(curvefilename, 'gdcurves')
  else
     deltaSIcurves=mpi_si2deltaSI_PD(sliceNum,seriesNum,outpath,numPreContrast_bld, numPreContrast_tiss, numSkip, scaleToStudyOne);
%     deltaSIcurves=deltaSIcurves(:,ranget); % added 10/6/05 so no need to run 1.1. again
     curvefilename=strcat(outpath,'deltaSIcurves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat');
     save(curvefilename, 'deltaSIcurves')
  end
drawnow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end;if ran(stg,4,'MODEL FITTING OF MRI DATA');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;if ran(stg,4.1,'Fit the regional MRI data to a 2-compartment model') ;
   modelType='full';
   fixedVp=99;  % cannot turn this off with full model
  % replaces bld curve with another one. What if scale factor different? Will subtract off and look the same... need to compare pre-contrast existing bld pool and new one? 
   %deltaSIcurves=mpi_useAlternateAIF(seriesNumAIF, sliceNumAIF,scaleAIF, sliceNum, seriesNum,outpath,numPreContrast_bld, numPreContrast_tiss, numSkip);
   h = figure;
%   deltaSIcurves=mpi_useAlternateAIF_tail(seriesNumAIF, sliceNumAIF,scaleAIF, sliceNum, seriesNum,outpath,numPreContrast_bld, numPreContrast_tiss, numSkip);%%% using tail frames to scale here...NP 082609
% TOOK OFF TAIL PART - using scale from diluting, dual bolus, 7/2011. 
   timeStampFileAIF=timeStampFile;
    if ~exist('framesToSelect','var')  % trying this EVRD 10/2/11
            framesToSelect=ranget;
    end
   deltaSIcurves=mpi_useAlternateAIF(seriesNumAIF, sliceNumAIF,scaleAIF, sliceNum, seriesNum,outpath,numPreContrast_bld, numPreContrast_tiss, numSkip,framesToSelect, ranget,timeStampFile, timeStampFileAIF);%%% using tail frames to scale here...NP 082609
   close(h);
   drawnow
   mycolors = lines(numAzimuthalRegions);
   %figure(seriesNum*10 + 2*sliceNum+1),subplot(2,2,4),cla, hold on, plot(deltaSIcurves(1,:),'k');
   figure(seriesNum*10 + 2*sliceNum+1),clf, hold on, plot(deltaSIcurves(1,:),'k');
   for region=1:numAzimuthalRegions
       plot(deltaSIcurves(region+1,:)','Color',mycolors(region,:));
       %pause
   end
% save with new blood
   if(lastFrame > size(deltaSIcurves,2) || lastFrame==0 || ~exist('lastFrame'))
       lastFrame = size(deltaSIcurves,2);
   end
   deltaSIcurves=deltaSIcurves(:,1:lastFrame);
 
   outfilename=strcat(outpath,'deltaSIcurves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.mat');
   save(outfilename, 'deltaSIcurves')
   %% run stage 4.2 right away:
   curvefilename=strcat(outpath,'deltaSIcurves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.mat');
   load(curvefilename);

   if exist(timeStampFile)
        fit_model_2comp_new(seriesNumAIF, sliceNumAIF, scaleAIF, sliceNum,seriesNum,outpath,'deltaSIcurves',delta_t,numAzimuthalRegions, numRadialRegions, flagPixelwise, 0, 0, fixedVp, fixedDelay, modelType)
   else
      %  fit_model_2comp(seriesNumAIF, sliceNumAIF, scaleAIF, sliceNum,seriesNum,outpath,'deltaSIcurves',delta_t,numAzimuthalRegions, numRadialRegions, flagPixelwise, 0, 0, fixedVp, fixedDelay, modelType, framesToSelect,ranget);
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;if ran(stg,4.2,'Show the 2-compartment model fits from stage 4.1');
  clusterFlag=0; modelType='full' %#ok<NOPRT>
     %infilename=strcat(outpath,'cinemri1.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat');
     %load(infilename);
     %mpi_displayRegionValues(cinemri1,sliceNum,seriesNum,outpath,numAzimuthalRegions, numRadialRegions, referenceFrame,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF, modelType)
     curvefilename=strcat(outpath,'deltaSIcurves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.mat');
     load(curvefilename);
     fixedVp=99;
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%% This is added to display the pixelwise maps of the values
     if(flagPixelwise);
        infilename=strcat(outpath,'cinemri1.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat');
        load(infilename);
        mpi_displayPixelwise(cinemri1,sliceNum,seriesNum,outpath,deltaSIcurves, referenceFrame, numPreContrast_bld, numPreContrast_tiss, useIntegralLinearFit,fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF,modelType)
       
     else
%      curvefilename=strcat(outpath,'deltaSIcurves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat');
%      load(curvefilename);   %EVRD 3/09, commented out, longer length than
%      one with .AIF in it... so lastFrame not being used... 

     infilename=strcat(outpath,'cinemri1.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat');
     load(infilename);
     
     if exist(timeStampFile)
         flagTimeStamps=1;
     else
         flagTimeStamps=0;
     end
     
%      if seriesNumAIF==-1 & flagPixelwise==0
%         mpi_displayRegionValues(cinemri1,sliceNum,seriesNum,outpath,numAzimuthalRegions, numRadialRegions, referenceFrame,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp)
%      else
         if flagPixelwise==0
        mpi_displayRegionValues(cinemri1,sliceNum,seriesNum,outpath,numAzimuthalRegions, numRadialRegions, referenceFrame,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF, modelType, flagTimeStamps)
     end
     end
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     showfits(deltaSIcurves,seriesNum,sliceNum,outpath,numAzimuthalRegions, numRadialRegions, 0,flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF, modelType, lastFrame, timeStampFile,timeStampFileAIF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;if ran(stg,4.3,'Compute upslope model') ;
    load(['Output/deltaSIcurves.study' int2str(seriesNum) '.slice' int2str(sliceNum) '.mat']);
    AIF = deltaSIcurves(1,:);
    st = length(AIF);
    upslopeRatios = zeros(size(deltaSIcurves,1),1);
    AIFmMax = 0;
    delays = 0;
    figure(seriesNum*10 + 2*sliceNum+1)
    subplot(2,2,4),clf,hold on,  
    slidingWindowWidth = 3; % Schwitter2001 used 3 point linear fit in blood and 5 point in tissue
    for t=1:(st-slidingWindowWidth)
        s = robustfit(((0:slidingWindowWidth)),AIF(t:(t+slidingWindowWidth)));
        if(s(2) > AIFmMax)
            AIFmMax = max(AIFmMax,s(2));
            upslopeRatios(1) = max(upslopeRatios(1),s(2));
            delays(1) = t;
        end
    end
    plot(AIF,'r','Linewidth',2.5);
    Y = AIF(delays(1):(delays(1)+slidingWindowWidth));
    X = delays(1):(delays(1)+slidingWindowWidth);
    %mpi_twoplot(Y,X,delays(1));
    slidingWindowWidth = 5;
    for curve=2:size(deltaSIcurves,1)
        for t=1:(st-slidingWindowWidth)
            s = robustfit(((0:slidingWindowWidth)),deltaSIcurves(curve,t:(t+slidingWindowWidth)));
            if(s(2) > upslopeRatios(curve))
                upslopeRatios(curve) = max(upslopeRatios(curve),s(2));
                delays(curve+1) = t;
            end
        end
        plot(deltaSIcurves(curve,:),'r','Linewidth',2.5);
        Y = deltaSIcurves(curve,delays(curve+1):(delays(curve+1)+slidingWindowWidth));
        X = delays(curve+1):(delays(curve+1)+slidingWindowWidth);
        %mpi_twoplot(Y,X,delays(curve+1));
    end
    upslopeRatios.m = upslopeRatios / AIFmMax;
    upslopeRatios.t = delays;
    save(['Output/upslopes.study' int2str(seriesNum) '.slice' int2str(sliceNum) '.mat'],'upslopeRatios');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  next from EVRD, in Sept14 mpi2d, commenting out other 4.31 below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;if ran(stg,4.31,'Compute upslope');
  
try   deltaSIcurves=mpi_useAlternateAIF(seriesNumAIF, sliceNumAIF,scaleAIF, sliceNum, seriesNum,outpath,numPreContrast_bld, numPreContrast_tiss, numSkip);
catch
    return;
end


% save with new blood
   orig_deltaSIcurves=deltaSIcurves;
   deltaSIcurves=deltaSIcurves(:,1:lastFrame);
  
   outfilename=strcat(outpath,'deltaSIcurves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.mat');
   save(outfilename, 'deltaSIcurves')

   curvefilename=strcat(outpath,'deltaSIcurves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.mat')
   load(curvefilename);


%   if  useDeltaSI~=1    %  12/30/03.
%      upslopes=mpi_upslope(sliceNum,seriesNum,outpath,'gdcurves', numAzimuthalRegions, numRadialRegions, flagPixelwise, numSkipUpslope);
%   else
     if exist(timeStampFile)
        timeStampFile
        upslopes=mpi_upslope(seriesNumAIF, sliceNumAIF, scaleAIF, sliceNum,seriesNum,outpath,'deltaSIcurves',delta_t, numAzimuthalRegions, numRadialRegions, flagPixelwise,numSkipUpslope, timeStampFile);
     else
         disp('no timeStampFile exists ')
        upslopes=mpi_upslope(seriesNumAIF, sliceNumAIF, scaleAIF, sliceNum,seriesNum,outpath,'deltaSIcurves',delta_t, numAzimuthalRegions, numRadialRegions, flagPixelwise,numSkipUpslope);
     end
%   end
  % should we use numPreContrast_bld and numSkip and use SI curves or [Gd]???

   deltaSIcurves=orig_deltaSIcurves;  % back to full length
   save(outfilename, 'deltaSIcurves')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % end;if ran(stg,4.31,'display Computed upslope model') ;
% %     
% %     infilename=strcat(outpath,'cinemri1.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat');
% %     load(infilename);
% %     load(['Output/upslopes.study' int2str(seriesNum) '.slice' int2str(sliceNum) '.mat']);
% %     AIF = upslopeRatios.m(1); %#ok<NASGU>
% %     tissue = upslopeRatios.m(2:end);
% %     delays = upslopeRatios.t; %#ok<NASGU>
% %     blood = load([outpath 'blood_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)]); %#ok<NASGU>
% %     endo = load([outpath 'endo_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)]);
% %     epi = load([outpath 'epi_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)]);
% %     angle = load([outpath 'Roi_start_angle.study' int2str(seriesNum) '.slice' int2str(sliceNum)]);
% % 
% % 
% %     %use the contours to recreate a mask as that's what the user would be
% %     %doing to recreate our results
% %     outtermask = roipoly(cinemri1(:,:,22),epi(:,1),epi(:,2));
% %     innermask = roipoly(cinemri1(:,:,22),endo(:,1),endo(:,2));
% %     newMyo = double((outtermask - innermask)>0);
% %     if(flagPixelwise)
% %         indicies = find(newMyo(:));
% %         for i=1:length(indicies)
% %             newMyo(i) = indicies(i);
% %         end
% %     else
% %         %use the myocardial mask to find the curves then the deltaSI curves
% %         [r,c] = find(outtermask>0);
% %         middle = [mean(r) mean(c)];
% %         [r,c] = find(newMyo>0);
% %         [theta,~] = cart2pol(r-middle(1), c - middle(2));
% %         theta = mod(theta - angle,2*pi);
% %         labeled = floor(theta/(2*pi)*numAzimuthalRegions);
% %         for i=1:length(r)
% %             newMyo(r(i),c(i)) = labeled(i)+1;
% %         end
% %     end
% %     [r,c] = find(newMyo>0);
% %     for i=1:length(r)
% %         newMyo(r(i),c(i)) = tissue(newMyo(r(i),c(i)));
% %     end
% %     figure(seriesNum*10 + 2*sliceNum+1), subplot(2,2,2), imagesc(newMyo), colorbar;title('Upslope of tissue/upslope of AIF');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;if ran(stg,5.1,'Fit the regional MRI data to a full Fermi model') ;
  % replaces bld curve with another one. What if scale factor different? Will subtract off and look the same... need to compare pre-contrast existing bld pool and new one? 
   modelType='fermiFull' %#ok<NOPRT>
   fixedVp=99;
   %deltaSIcurves=mpi_useAlternateAIF(seriesNumAIF, sliceNumAIF,scaleAIF, sliceNum, seriesNum,outpath,numPreContrast_bld, numPreContrast_tiss, numSkip);
   deltaSIcurves=mpi_useAlternateAIF_tail(seriesNumAIF, sliceNumAIF,scaleAIF, sliceNum, seriesNum,outpath,numPreContrast_bld, numPreContrast_tiss, numSkip);%%% using tail frames to scale here...NP 082609
   drawnow

% save with new blood
   deltaSIcurves=deltaSIcurves(:,1:lastFrame);
 
   outfilename=strcat(outpath,'deltaSIcurves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.mat');
   save(outfilename, 'deltaSIcurves')
   %% run stage 4.2 right away:
%   curvefilename=strcat(outpath,'deltaSIcurves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.mat')
   curvefilename=strcat(outpath,'deltaSIcurves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.mat');
   load(curvefilename);

   if exist(timeStampFile)
        fit_modelFermi(seriesNumAIF, sliceNumAIF, scaleAIF, sliceNum,seriesNum,outpath,'deltaSIcurves',delta_t,numAzimuthalRegions, numRadialRegions, flagPixelwise, 0, 0, fixedVp, fixedDelay, modelType, timeStampFile, timeStampFileAIF);
   else
        fit_modelFermi(seriesNumAIF, sliceNumAIF, scaleAIF, sliceNum,seriesNum,outpath,'deltaSIcurves',delta_t,numAzimuthalRegions, numRadialRegions, flagPixelwise, 0, 0, fixedVp, fixedDelay, modelType);
   end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;if ran(stg,5.2,'Show the Fermi model fits from stage 5.1');
  clusterFlag=0; modelType='fermiFull' %#ok<NASGU,NOPRT>
  fixedVp=99;
  %   infilename=strcat(outpath,'cinemri1.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat');
  %   load(infilename);
  %   mpi_displayRegionValues(cinemri1,sliceNum,seriesNum,outpath,numAzimuthalRegions, numRadialRegions, referenceFrame,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF, modelType)
  %   curvefilename=strcat(outpath,'deltaSIcurves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.mat')
     curvefilename=strcat(outpath,'deltaSIcurves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.mat');
     load(curvefilename);
     showfits(deltaSIcurves,seriesNum,sliceNum,outpath,numAzimuthalRegions, numRadialRegions, 0,flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF, modelType, lastFrame, timeStampFile,timeStampFileAIF);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;if ran(stg,6.1,'Show the k-trans values in a pretty plot');
    template = ['.study' num2str(seriesNum) '.slice' num2str(sliceNum)];
    infilename=strcat(outpath,'cinemri1.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat');
    load(infilename);
    [sx,sy,~] = size(cinemri1);
    flowFiles = dir(['Output/flowvalues' template '.6.1_fixedDelay0.AIF_' num2str(seriesNum) '_' num2str(sliceNum) '_1.txt.full.txt']);
    flowFiles = dir(['Output/flowvalues' template '.6.1_fixedDelay0.AIF_' num2str(seriesNumAIF) '_' num2str(sliceNumAIF) '_',num2str(scaleAIF),'.txt.full.txt']);
    if(isempty(flowFiles))
        disp('empty flowFiles, checkin tino this ')
        keyboard;
    else
        clear ktrans ve
        ktransValues = zeros(numAzimuthalRegions,1);   %should grow if more regions are encountered
        ve = zeros(numAzimuthalRegions,1);
        fid = fopen(['Output/' flowFiles(1).name]);
        tline = fgetl(fid);
        while(ischar(tline))
            if(~isempty(strfind(tline,'Ktrans')))
                A = sscanf(tline,'Ktrans  %d     %f');
                ktransValues(A(1)) = A(2);
            end
            if(~isempty(strfind(tline,'ve')))
                A = sscanf(tline,'ve    %d      %f');
                ve(A(1)) = A(2);
            end
            tline = fgetl(fid);
        end
        fclose(fid);
        
        ktransmin = min(ktransValues); %#ok<NASGU>
        ktransmax = max(ktransValues); %#ok<NASGU>
        myscaledktrans = ktransValues;
        myscaledktrans = myscaledktrans / max(myscaledktrans);
        
        %construct radialWeight
        blood = load([outpath 'blood_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)]); %#ok<NASGU>
        endo = load([outpath 'endo_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)]);
        epi = load([outpath 'epi_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)]);
        angle = load([outpath 'Roi_start_angle.study' int2str(seriesNum) '.slice' int2str(sliceNum)]);

        legacy = 0;
        if(angle > 2*pi)
            legacy = 1;
            angle = mod((angle-90)*pi/180,2*pi);
        end
        
        %use the contours to recreate a mask as that's what the user would be
        %doing to recreate our results
        outtermask = roipoly(cinemri1(:,:,22),epi(:,1),epi(:,2));
        [r,c] = find(outtermask>0);
        middle = [mean(r) mean(c)];
        innermask = roipoly(cinemri1(:,:,22),endo(:,1),endo(:,2));
        newMyo = (outtermask - innermask)>0;
        [r,c] = find(newMyo>0);
        [theta,radius] = cart2pol(r-middle(1), c-middle(2));
        theta = mod(theta - angle,2*pi);
        
        myregionalMask = zeros(sx,sy);
        labeled = floor(theta/(2*pi)*length(ktransValues));
        %make it so that  region 3 and 4 are in the septal wall
        %(the region counter clockwise from the insertion angle is region 3, and the next clockwise region is region 4)
        if(~legacy)
            labeled = mod(labeled + 2,6);
        end
        for i=1:length(r)
            myregionalMask(r(i),c(i)) = labeled(i)+1;
        end
        
        radialWeight = zeros(sx,sy);
        dtheta = .07;
        for th=0:.1:(2*pi)

            indicies = (theta <= (th+dtheta)).* (theta >= (th-dtheta));
            meanRadius = mean(radius(indicies>0)); %#ok<%#ok<MSNU> NASGU>
            stdRadius = std(radius(indicies>0))*MyoFatness; %#ok<%#ok<MSNU> NASGU>
            myi = find(indicies>0);
            for i=1:length(myi)
                radialWeight(r(myi(i)),c(myi(i))) = Alpha;
                radialWeight(r(myi(i)),c(myi(i))) = Alpha*exp(-(radius(myi(i)) - meanRadius)^2/(2*stdRadius^2));
            end
        end
        %radialWeight = filter2(fspecial('gaussian',[sx,sy],1),radialWeight);
    
        %temp = mean(cinemri1(:,:,(22-3):(22+3)),3);
        temp = cinemri1(:,:,referenceFrame);
        temp = temp - min(temp(:));
        temp = temp / max(temp(:));
        myBW = temp;
        myrgb = zeros(sx,sy,3);
        myrgb(:,:,1) = temp;myrgb(:,:,2) = temp;myrgb(:,:,3) = temp;
        load('spect.cmap');
        mycolors = spect/max(spect(:)); %#ok<NODEF>
        ncolors = length(mycolors);
        mycolors = reshape(mycolors,1,ncolors,3);
        for x=1:sx
            for y=1:sy
                if(myregionalMask(x,y) > 0)
                    c1 = mycolors(1,1+round((ncolors-1)*myscaledktrans(myregionalMask(x,y))),:);
                    c2 = myrgb(x,y,:);
                    myrgb(x,y,:) = radialWeight(x,y)*c1 + (1-radialWeight(x,y))*c2;
                end
            end
        end
        
        %display the k-transes in a seperate figure
        figure(seriesNum); clf
        set(gcf,'position',[50 600 600 300])
        if(FinalShading)
            [X,Y] = meshgrid(1:sx,1:sy);
            tcolor(X,Y,myrgb);
            axis ij
        else
            imagesc(myrgb)
        end
        hold on,colormap(spect/max(spect(:))), cb = colorbar;
        temp = infile(1:(end-1));
        i = strfind(temp,'/');
        set(gcf,'Name',temp((i(end)+1):end));
        [sortedScaledktrans IX] = sort(myscaledktrans); %#ok<ASGLU>
        sortedktrans = ktransValues(IX);
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        %set(cb,'YTick',[0; sortedScaledktrans]);
        set(cb,'YTickLabel',{round(10*(0:max(sortedktrans)/10:max(sortedktrans)))/10});
        myvalues = myregionalMask(myregionalMask > 0);
        [x,y] = find(myregionalMask > 0);
        x = x - middle(1); y = y - middle(2);
        [theta,radius] = cart2pol(x,y);
        theta = mod(theta,2*pi);
        for region=1:length(ktransValues)
            thetas = theta(myvalues == region);
            radii = radius(myvalues == region);
            hadToWrap = 0;
            if(any(thetas > (2*pi - pi/4)) && any(thetas < pi/4))
                thetas = mod(thetas + pi,2*pi);
                hadToWrap = 1;
            end
            mytheta = median(thetas);
            myr = 1.3*max(radius(myvalues == region));
            if(hadToWrap)
                [myx,myy] = pol2cart(mytheta-pi,myr);
            else
                [myx,myy] = pol2cart(mytheta,myr);
            end
            if(~exist('HorzOffset'))
                HorzOffset = 0;
            end
            myx = myx + middle(1);myy = myy + middle(2);
            %test target location for color
            tx = max(2,min(round(myx),sx-1)); tx = (tx-1):(tx+1);
            ty = max(2,min(round(myy),sy-1)); ty = (ty-1):(ty+1);
            
            mybackColor = mean(mean(myBW(tx,ty)));
            mybackColor = mybackColor<.5;
            
            text(myy+HorzOffset,myx,['Reg. ' num2str(region) ' = ' num2str(round(10*ktransValues(region))/10)],'FontSize',12,'Color',[mybackColor mybackColor mybackColor],'FontName','Helvetica');
            rightHandSideThetas = thetas(thetas < (min(thetas)+pi/30));
            rmin = min(radii(thetas < (min(thetas)+pi/30)));
            rmax= max(radii(thetas < (min(thetas)+pi/30)));
            if(hadToWrap)
                [x,y] = pol2cart([min(rightHandSideThetas) min(rightHandSideThetas)]-pi,[rmin rmax]);
            else
                [x,y] = pol2cart([min(rightHandSideThetas) min(rightHandSideThetas)],[rmin rmax]);
            end
            x = x + middle(1);y = y + middle(2); %#ok<NASGU>
            
            %line(y,x,'Color',[1 .5 0]);
        end
        extendedEndo = vertcat(endo((end-9):end,:),endo,endo(1:10,:));
        extendedEpi = vertcat(epi((end-9):end,:),epi,epi(1:10,:));
        clear smoothedExtendedEndo smoothedExtendedEpi
        smoothedExtendedEndo(:,1) = conv(fspecial('gaussian',[length(extendedEndo) 1],.5),extendedEndo(:,1),'same');
        smoothedExtendedEndo(:,2) = conv(fspecial('gaussian',[length(extendedEndo) 1],.5),extendedEndo(:,2),'same');
        smoothedExtendedEpi(:,1) = conv(fspecial('gaussian',[length(extendedEpi) 1],.5),extendedEpi(:,1),'same');
        smoothedExtendedEpi(:,2) = conv(fspecial('gaussian',[length(extendedEpi) 1],.5),extendedEpi(:,2),'same');
        smoothedEndo = smoothedExtendedEndo(11:(end-10),:); %#ok<NASGU>
        smoothedEpi = smoothedExtendedEpi(11:(end-10),:); %#ok<NASGU>
        plot(smoothedEndo(:,1), smoothedEndo(:,2));
        plot(smoothedEpi(:,1), smoothedEpi(:,2));
        title('Regional K-trans');
        drawnow
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(wordy), disp(['Time=' num2str(toc) ' seconds']); end
return;














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   COMMON SERVICE ROUTINES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert arguments 
% into a character buffer
% 

function buf=arg2buf(n,v)

buf=[];
for k=1:n
    if(ischar(v{k}))
        buf=[buf v{k} sprintf('\n')];
    end
end

% convert character buffer
% to a structure array
%
function par=buf2par(par,buf)

if(isempty(par)), return; end;
lf=sprintf('\n');
ind=[];
if(isempty(buf)==0)
 ind=find((buf==lf)|(buf==';'));
end

if(isempty(ind)==0), buf(ind)=' '; end;

ind=[1 ind length(buf)];
ii=size(par,1);

for k=1:length(ind)-1
   str=buf(ind(k):ind(k+1));
   in=find(str=='%', 1 );
   if(~isempty(in)), str=str(1:in-1); end;
   in=find(str=='=');
   if(~isempty(in))
      ii=ii+1;
      par{ii,1}=delspace(str(1:in-1));
      par{ii,2}=str(in+1:end);
   end
end

return

% function res=delspace(str)

% in=find(str~=' ');
% if(length(in)>0) res=str(in); end;

% delete spaces from the string

function res=delspace(str)

 in=find(str~=' ');
 if(~isempty(in)), res=str(in); end;

% Grabbing the parameter 
% from the input buffer
%
%
function res=getpar(pars,kind,name,def,cmnt)

if(isempty(pars))
   fprintf('%8s = %-12s  %s',name,def,cmnt)
   res=[]; 
else
   in=find(strcmp(pars(:,1),name), 1 );
   res=char(pars(in,2));
end

if(isempty(res)), res=def; end;
if(isempty(res)),  return; end;

switch(lower(kind))
  case 'char' 
    if(ischar(res)==0), res=num2str(res); end; 
    res=strtrim(res);
  case 'double'
    res=txt2double(res);
end;

return;


% run - interface function
%       
%    flag=ran(stg,stgnum,msg)
%
%
%
function flag=ran(stg,stgnum,msg)

in=find(abs(stg(3:end)-stgnum)<0.00001, 1);

if(~isempty(in))            % stage touched
 if(stg(1)), disp([num2str(stgnum) ' ' msg]); end; % wordy
 flag=stg(2);
else
 flag=0;
end

% READTXT read file into char buffer
% if no file, return empty
% 
%  buf= readtxt(in)
%
function buf=readtxt(in)
buf=[];

ff=fopen(in,'r');
if(ff==-1), return; end;

while 1
  line = fgetl(ff);
  if ~ischar(line), break, end
  buf=[buf line sprintf('\n')];
end
fclose(ff);

return;

function deliver(codename,dir1)

codename=[codename '_all.m'];
lst1=dir([dir1 '/*.m']);

cmd=['function t=codedate; t=''' date ''';'];
ii=fopen([dir1 '/codedate.m'],'w');
fprintf(ii,'%s\n',cmd);
fclose(ii);

k=1;
lst{k}=[dir1 '/' codename];

for l=1:length(lst1(:))
 if(length(findstr(lst1(l).name,codename ))==0) %#ok<REMFF1,ISMT>
  k=k+1;
  lst{k,1}=[ dir1 '/' lst1(l).name ];
 end
end

dos(['cat ' lst{1} ' >' codename ]);
for l=2:k
 cmd=['cat ' lst{l} ' >>' codename ];
 dos(cmd);
end

%disp('converting to pcode, can take a  minute, or often hangs!!!')
%pcode(codename);
%dos(['rm ' codename]);
%disp('done converting to pcode')

return

% dos('rm *.p');
% pcode(codename);
% dos(['rm '  codename]);
%  PS2LATEX Convert matlab postsript file to latex-compartable format
%
%           ps2latex(psfile);
%
%           Reads postscript file bounding box from the end
%           of file and puts it to the file beginning, 
%           so that this file can be incorporated to latex 
%           documents as a figure.
%           
%  Example: 
%           plot([1 2 3]);
%           print -dps fig.ps
%           ps2latex('fig.ps');
%
% See also: PRINT
%

% to speed up, swallow the whole file, change and output for new
% close
% clear
% id=fopen('tst.m','r');
% p=fscanf(id,'%c');
% fclose(id);
% p


%
%  02.27.97
%  (c) 1997 Oleg Portniaguine
%
%
function ps2latex(psfile) %#ok<DEFNU>

% open and position near the end

inp=fopen(psfile,'r');
fseek(inp,-100,'eof');

% find pages and bounding box values

while(1)
%  [s,n]=readstr(inp);
%  if n==0; break; end;
 
 s=fgetl(inp);
 if ~ischar(s), break, end;

 if( isempty(strfind(s,'%%Pages:')) ==0)
     pg=s(9:length(s));
 elseif( isempty(strfind(s,'%%BoundingBox:')) ==0)
     bb=s(15:length(s));
 end
end

if (isempty(pg)==1) 
  error('Cannot find number of pages at the end of file');
end
if (isempty(bb)==1) 
  error('Cannot find bounding box at the end of file');
end

% read and change the beginning of file

frewind(inp);
ss=[];

while(isempty(bb)==0)

%[s,n]=readstr(inp);
%if n==0; break; end;

 s=fgetl(inp);
 if ~ischar(s), break, end;

 if( strcmp(s,'%%Pages: (atend)') ==1)
   s=['%%Pages: ' pg];
 elseif( strcmp(s,'%%BoundingBox: (atend)') ==1)
   s=['%%BoundingBox:' bb];
   bb=[];
 end

 ss=[ss sprintf('%s\n',s)];

end

% read the rest and close

s=fread(inp);
fclose(inp);

% rewrite the file

out=fopen(psfile,'w');
fwrite(out,abs(ss),'uchar');
fwrite(out,s,'uchar');
fclose(out);
%
% Initial help guide for the code
%
function dispHelp(code)

 disp(' ');
 disp('To list all parameters and their default values:');
 disp([ ' ' code ' lst=par']);
 disp('To list all stages from 1 to 10:');
 disp([ ' ' code ' lst=all']);
 disp('To run stage 1.1  :');
 disp([ ' ' code ' stg=1.1 ']);
 disp(' ');

return


function cinemri1 = RetrieveMovieFromMat(infile,rangex, rangey,ranget,flipy,flipx,rotNinety,inendi,intype,inpnxy)
      if (infile(length(infile)-3:end) ~= '.mat')
          ff=fopen(infile,'r',inendi);
          a=fread(ff,intype);
          fclose(ff);
          n3=length(a)/prod(inpnxy);
          a=double(reshape(a,[inpnxy n3]));
          if(flipx)
              cinemri=a(max(rangex):-1:min(rangex),rangey,ranget); % this drops out unwanted frames but need to worry about delta_t!!!!
          else
              cinemri=a(rangex,rangey,ranget);
          end
          if(rotNinety)
              for ii=1:max(ranget)-min(ranget)+1
                  cinemri(:,:,ii)=rot90(cinemri(:,:,ii),-1);
              end
          end

          for i=1:(max(ranget)-min(ranget)+1)         % transpose for display purposes
              tmpp(:,:,i)=cinemri(:,:,i)';
          end
          cinemri1=tmpp;

      else
          filename=infile; %%% assuming the images are in .mat format and can be loaded directly
          load (filename)
          a=mpi_upsampleImages(imgrr(rangey,rangex,ranget));

          if(flipy)
              cinemriy=a(end:-1:1,:,:); % this drops out unwanted frames but need to worry about delta_t!!!!
          else
              cinemriy=a(:,:,:);
          end
          %%% often the reconstructed .mat images require 'flipping' in the x- and y-directions
          if(flipx)
              cinemri3=cinemriy(:,end:-1:1,:); % this drops out unwanted frames but need to worry about delta_t!!!!
          else
              cinemri3=cinemriy;
          end


          if(rotNinety)
              for ii=1:max(ranget)-min(ranget)+1
                  cinemri(:,:,ii)=rot90(cinemri3(:,:,ii),-1);
              end
          else
              cinemri(:,:,:)=cinemri3(:,:,:);
          end

          for i=1:(max(ranget)-min(ranget)+1)         % transpose for display purposes
              tmpp(:,:,i)=cinemri(:,:,i);
          end
          cinemri1=tmpp;
      end
% conversion from string to double
% 
% colon conversion here
% vectors without colons are converted anyway by str2num
% 
function res=txt2double(res)
            
if(ischar(res)) 
       
       in=find(res==':');
       inb=find( (res=='[')+(res==']') );

       switch(length(in))
        case 0
         res=str2num(res);
         return;
        case 1
         res(inb)=' ';         
         stp=1;
         beg=str2num(res(1:in-1));
         fin=str2num(res(in+1:end));
        case 2
         res(inb)=' ';
         beg=str2num(res(1:in(1)-1));
         stp=str2num(res((in(1)+1):(in(2)-1)));
         fin=str2num(res((in(2)+1):end));
        otherwise
        disp(['Error in parameter ' name ' bad value ' res ]); 
       end
       res=[ beg : stp : fin ]; %#ok<NBRAK>
end
return;
       




%
% Compute static shift for 
% the whole image
% grids=[4 2 1]
%
function [m2,shfts]=staticshft(wordy,m,grids,referenceFrame, numSkipRegistration) %#ok<DEFNU>

%if ~exist(numSkipRegistration)
%   numSkipRegistration=0;
%end
% new weights  from Ganesh, 8/08/05:
ww=0.5+cosw2d(size(m)); % weighting function
%cc=chebwin(100,60);
%cc=chebwin(100,50);
%ww=cc*cc';
[nx,ny]=size(m); %#ok<NASGU,ASGLU>
%ww=1;  % didn't seem to make a difference for the 1 dataset I looked at carefully - ED 11/4/03  
m2=m;
w3=zeros(5,5);

 mframe=m2(:,:,referenceFrame);  % do not update reference frame after
% mframe=mean(m2,3);  % do not update reference frame after
                     % each iteration

for it=1+numSkipRegistration:size(m2,3)


    sh=[0 0];
    for n=grids

        [x,y]=ndgrid(sh(1)+n*[-2:2],sh(2)+n*[-2:2]); %#ok<NBRAK>

        for k=1:length(x(:))
            r=imgshft(m2(:,:,it),[x(k) y(k)])- mframe;
            %    r = shiftimage(m2(:,:,it),y(k)/ny,x(k)/nx) - mframe ;
            r=r.*ww;
            w3(k)=r(:)'*r(:);
        end

        [~,i1]=min(w3(:));
        sh=[x(i1) y(i1)];

        %debugging plots
        %clf;
        %contourf(x,y,w3); hold on; plot(sh(1),sh(2),'r*');
        %title('Misfit');
        %colorbar;
        %pause;
    end

    if(wordy>=0), disp([it sh]); end
    m2(:,:,it)=imgshft(m2(:,:,it),sh);
    % m2(:,:,it) =shiftimage(m2(:,:,it),sh(2)/ny,sh(1)/nx) ;
    shfts(it,:)=sh;

end


% shift the image on sh pixels
%
function m2=imgshft(m,sh)

m2=intshft(m,floor(sh));
% don't use below unless have full FOV, and may want to modify so only wraps in
% phase encode direction
%m2=intshftCircular(m,floor(sh));
 sh=sh-floor(sh);
 if(sh==0), return; end;
% m2=remshft(m2,sh);
% Interpolation
return; 

%
function m2=remshft(m,sh) %#ok<DEFNU>

[nx,ny]=size(m);
m2=m;
m2(1:nx-1,1:ny-1)= ...
(m2(1:nx-1,1:ny-1)*(1-sh(1)) + m2(2:nx,1:ny-1)*sh(1)) * (1-sh(2)) + ...
(m2(1:nx-1,2:ny)*(1-sh(1)) + m2(2:nx,2:ny)*sh(1)) * sh(2); %#ok<NASGU>

m2=shiftimage(m,sh(2)/ny,sh(1)/nx) ;

%mfft=fft2(m);
%[M N]= size(mfft);

%mtmp=fftshift(mfft);
%mtmp=mfft;
%mfft=[mtmp zeros(M,3*N)];  mfft = [mfft; zeros(3*M,4*N)];
%mfft=[zeros(M,2*N) mtmp zeros(M,2*N)];  mfft = [zeros(2*M,5*N); mfft; zeros(2*M,5*N)];
%mfft=fftshift(mfft);
% moves to right sh(1) and up sh(2) , for integer 
%sh=5*sh
%shiftkx=exp(-i*2*pi*sh(2)*(1:5*M)/(5*M)); shiftky=exp(-i*2*pi*sh(1)*(1:5*N)/(5*N));
%xcenter=(M-1)/2;
%shiftkx=exp(-i*2*pi*sh(2)*(0:(M/2-1))/(M)); shiftkx=[shiftkx exp(-i*2*pi*sh(2)*(((M/2-1):-1:0)/M))];
%shiftky=exp(-i*2*pi*sh(1)*(0:(N/2-1))/(N)); shiftky=[shiftky exp(-i*2*pi*sh(1)*(((N/2-1):-1:0)/N))];


%shiftkx=fftshift(exp(-i*2*pi*sh(2)*(0:(M-1))/(M)));
%shiftky=fftshift(exp(-i*2*pi*sh(1)*(0:(N-1))/(N)));
%ll=fft(m,ny,2);
%tmpp=zeros(size(m));
%for n=1:ny
%   tmpp(:,n)=ll(:,n).*shiftkx';
%end
%m2=ifft(tmpp,ny,2);
%ll=fft(m2,nx,1);
%tmpp=zeros(size(m));
%for n=1:nx
%   tmpp(n,:)=ll(n,:).*shiftky;
%end
%m2=ifft(tmpp,nx,1);

%%shiftsfourier=fftshift(shiftkx'*shiftky);
%%shiftsfourier=(shiftkx'*shiftky);
%%mfft = mfft.*shiftsfourier;
%%m2=ifft2(mfft,M,N);
%imagesc(abs(m2));

return;

% Purely integer shift
%
function m2=intshft(m,sh)

[nx,ny]=size(m);

inx=(1+abs(sh(1))) : (nx-abs(sh(1)));
iny=(1+abs(sh(2))) : (ny-abs(sh(2)));

m2=m;
m2(inx-sh(1),:)= m(inx,:);      % this is original
m2(:,iny-sh(2))= m2(:,iny);

%disp('in intshift')
%min(min(m))
%min(min(m2))

%m2=zeros(size(m));
%m2(inx-sh(1),:)= m(inx,:);      % could make this different?
%m3=zeros(size(m));
%m3(:,iny-sh(2))= m2(:,iny);
%m2=m3;

return;

function imgnew=intshftCircular(img,sh) %#ok<DEFNU>
% above function cuts off both ends (just leaves as existing m).^M
% add 3/02 by ED to wrap ( circular shifts - will need for unfold)^M
% may want to only wrap in PE direction?

        xshift=-sh(1); yshift=-sh(2);
        [xdim, ydim]=size(img);
        imgnew=img;
        if xshift > 0
                imgnew(1+xshift:xdim,:)=img(1:xdim-xshift,:);
                imgnew(1:xshift,:)=img(xdim-xshift+1:xdim,:);
        end
        if xshift < 0
                imgnew(1:xdim+xshift,:)=img(1-xshift:xdim,:);
                imgnew(xdim+xshift+1:xdim,:)=img(1:-xshift,:);
        end
   	  img=imgnew;
        if yshift > 0
                imgnew(:,1+yshift:ydim)=img(:,1:ydim-yshift);
                imgnew(:,1:yshift)=img(:,ydim-yshift+1:ydim);
        end
        if yshift < 0
                imgnew(:,1:ydim+yshift)=img(:,1-yshift:ydim);
                imgnew(:,ydim+yshift+1:ydim)=img(:,1:-yshift);
        end

return;
