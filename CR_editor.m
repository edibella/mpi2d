function varargout = CR_editor(varargin)
% CR_EDITOR M-file for CR_editor.fig
%      CR_EDITOR, by itself, creates a new CR_EDITOR or raises the existing
%      singleton*.
%
%      H = CR_EDITOR returns the handle to a new CR_EDITOR or the handle to
%      the existing singleton*.
%
%      CR_EDITOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CR_EDITOR.M with the given input arguments.
%
%      CR_EDITOR('Property','Value',...) creates a new CR_EDITOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CR_editor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CR_editor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CR_editor

% Last Modified by GUIDE v2.5 22-Dec-2010 09:54:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CR_editor_OpeningFcn, ...
                   'gui_OutputFcn',  @CR_editor_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before CR_editor is made visible.
function CR_editor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CR_editor (see VARARGIN)

% Choose default command line output for CR_editor
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% UIWAIT makes CR_editor wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global showImageMutex bumpingMutex dragMutex keyPressMutex v
if(exist('v','var') && isstruct(v))
    names = fieldnames(v);
    v = rmfield(v,names);
end
keyPressMutex = 0;
dragMutex = 0;
showImageMutex = 0;
bumpingMutex = 0;
%v = get(handles.figure1,'UserData');


v.bumping = 0;
v.r = 6;
v.playing = 0;
v.showSmooth = 0;
v.OldTimeStamp = now;
v.changedShifts = 0;
v.HitAKeyTimeStamp = now;
v.BackgroundImageType = 'Frames';
v.DeletePressed = 0;
v.changedContours = 0;
v.StartAngle = 3*pi/2;
v.blood = [];
v.legacy = 0;
v.waitingForAngleClick = 0;
v.mouse = [1 1];v.oldmouse = [1 1];
v.clickType = '';
v.showMyoMask = 0;
%% fix the axes so the markers are either missing or on the correct side
axes(handles.axes1); %#ok<MAXES>
set(gca,'XTick',[]);
set(gca,'YTick',[]);

axes(handles.axes2); %#ok<MAXES>
set(gca,'YAxisLocation','right');
set(gca,'YGrid','on');
set(gca,'YMinorGrid','on');
set(gca,'XMinorGrid','off');
axes(handles.axes3); %#ok<MAXES>
set(gca,'XTick',[]);
set(gca,'YAxisLocation','right');
axes(handles.axes4); %#ok<MAXES>
set(gca,'XTick',[]);
set(gca,'YAxisLocation','right');

%% parse the parameters
if(isempty(varargin))
    cd('/v/raid1/bmatthew/MRIdata/Cardiac/Brian-Testing/P050520A/Processing');
    varargin = {61000,1};
end
if(length(varargin)==2)
    if(varargin{1} == -1)
        disp('No valid mpi2d.par file')
        return;
    end
    v.seriesNum = varargin{1};
    
    v.sliceNum = varargin{2};
    set(handles.figure1,'Name',['CR_editor( Series: ' num2str(v.seriesNum) ' - Slice: ' num2str(v.sliceNum) ' )']);
    v.parfileName = ['series' num2str(v.seriesNum) '.slice' num2str(v.sliceNum) '.par'];
    ParFileName = v.parfileName; %#ok<NASGU>
    ReadPar;
    referenceFrame=17  %EVRD 4/27/11, HACK!  where does this come from?? tryied putting in .par
    v.ContourData.ReferenceFrame = referenceFrame;
    
    %set(handles.figure1,'UserData',v);
    sucessful = loadCR(handles);
    if(~sucessful)
        return;
    end
    %v = get(handles.figure1,'UserData');
    
    
    v.t = 1;
    for i=1:3
        v.contourSelected = i;
        %set(handles.figure1,'UserData',v);
        MasksChanged(handles);
        %v = get(handles.figure1,'UserData');
    end

    %set(handles.figure1,'UserData',v);
    set(handles.radiobutton3,'Value',v.contourSelected == 3);
    set(handles.radiobutton4,'Value',v.contourSelected == 2);
    set(handles.radiobutton5,'Value',v.contourSelected == 1);
    
    showImage(handles);
    showCurves(handles);
    showShifts(handles);

else
    disp('Please pass in the series and slice number');
    return;
end
function sucessful = loadCR(handles,UpdateOnlyRecent)
global v
    if(~isfield(handles,'figure1')), return, end
    if(~exist('UpdateOnlyRecent','var'))
        UpdateOnlyRecent = 0;
    end
    %v = get(handles.figure1,'UserData');
    sucessful = 0;
    %movie
    temp = dir(['Output/cinemri1.study' num2str(v.seriesNum) '.slice' num2str(v.sliceNum) '.mat']);
    if(isempty(temp))
        disp('Cine file was not created.  Have you not processed this dataset?');
        return;
    end
    if(~UpdateOnlyRecent || temp.datenum > v.OldTimeStamp)
        load(['Output/cinemri1.study' num2str(v.seriesNum) '.slice' num2str(v.sliceNum) '.mat']);
        
        %movie
        v.DCData.DCmri = cinemri1;
        [v.sx v.sy v.st] = size(v.DCData.DCmri);
        for t=1:v.st
            v.limits(t,:) = [min(min(v.DCData.DCmri(:,:,t))); max(max(v.DCData.DCmri(:,:,t)))];
        end
        v.globalLimits = [min(v.DCData.DCmri(:)) max(v.DCData.DCmri(:))];
        if(isfield(v,'delays'))
            v = rmfield(v,'delays');
        end
        if(isfield(v,'Manhattan'))
            v = rmfield(v,'Manhattan');
        end
        if(isfield(v,'Upslope'))
            v = rmfield(v,'Upslope');
        end
    end

    
    %curves
    temp = dir(['Output/curves.study' num2str(v.seriesNum) '.slice' num2str(v.sliceNum) '.mat']); 
    if(isempty(temp))
         v.CurvesData.NAzimuthalRegions = 6;
         % someday handle this gracefully...   this is incomplete now
         % 9/19/11 EVRD
    else
         
     if(~UpdateOnlyRecent || temp.datenum > v.OldTimeStamp)
        temp = load(['Output/curves.study' num2str(v.seriesNum) '.slice' num2str(v.sliceNum) '.mat']);
        
        v.CurvesData.NAzimuthalRegions = size(temp.curves,1)-1;
        for i=1:v.CurvesData.NAzimuthalRegions
            v.CurvesData.curves(i,:) = temp.curves(i+1,:);
        end
        v.AIFData.AIF = temp.curves(1,:);
     end
    end
    
    
    %DeltaSIcurves
    temp = dir(['Output/deltaSIcurves.study' num2str(v.seriesNum) '.slice' num2str(v.sliceNum) '.mat']);
    
    if(isempty(temp))
        % someday handle this gracefully...   v.blood = [];
    else
        
        if(~UpdateOnlyRecent || temp.datenum > v.OldTimeStamp)
        temp = load(['Output/deltaSIcurves.study' num2str(v.seriesNum) '.slice' num2str(v.sliceNum) '.mat']);
        v.blood = temp.deltaSIcurves(1,:);
        v.CurvesData.NAzimuthalRegions = size(temp.deltaSIcurves,1)-1;
        for i=1:v.CurvesData.NAzimuthalRegions
            v.CurvesData.deltaSIcurves(i,:) = temp.deltaSIcurves(i+1,:);
        end
        v.AIFData.deltaSIAIF = temp.deltaSIcurves(1,:);
        end
    
    
    
    v.KtransData = loadMostRecentKtrans(v.seriesNum,v.sliceNum);
    end
    
    %shifts
    temp = dir(['Output/shiftsMAN.study' num2str(v.seriesNum) '.slice' num2str(v.sliceNum) '.txt']);
    temp2 = dir(['timeStampSer' num2str(v.seriesNum) '.mat']);
    if(~UpdateOnlyRecent || temp.datenum > v.OldTimeStamp || temp2.datenum > v.OldTimeStamp)
        shifts = load(['Output/shiftsMAN.study' num2str(v.seriesNum) '.slice' num2str(v.sliceNum) '.txt']);
%         tmp_label=v.seriesNum - 1000*(floor(v.seriesNum/1000)); %EVRD
%         if tmp_label >= 90 && tmp_label <= 96
%             load(['timeStampSer' num2str(v.seriesNum - tmp_label) '.mat']);  %round
%         else
            %load(['timeStampSer' num2str(v.seriesNum) '.mat']); 
        %end
        try    
            load(['timeStampSer' num2str(v.seriesNum) '.mat']);          
        catch
            load(['timeStampSer' num2str(1000*(floor(v.seriesNum/1000))) '.mat']);  %round  EVRD
        end
        
        if(length(shifts) ~= size(v.DCData.DCmri,3))
            disp('Length of shift vector does not equal the duration of the movie');
            disp('I leave it in your capable hands to resolve it and save the data correctly so when it''s loaded everything is consistant');
            keyboard;
        end
        if(length(timeStamp) ~= size(v.DCData.DCmri,3))
            disp('Length of timestamp vector does not equal the duration of the movie');
            disp('Now chopping first entry of timestamp, since assuming recons skipped first frame. Need to handle this in stage 4.1 too! ');
            timeStamp=timeStamp(2:end);
            %keyboard;
        end
        v.timeStamp = timeStamp;
        
        for t=1:size(v.DCData.DCmri,3)
            v.ShiftsData(t).shift = shifts(t,:);
            v.ShiftsData(t).timestamp = v.timeStamp(t);
            v.DCData.timestamp(t) = v.timeStamp(t);
        end
    end

    %contours
    temp = dir(['Output/blood_polyCoords.study' num2str(v.seriesNum) '.slice' num2str(v.sliceNum)]);
    if(isempty(temp))
        %obviously the contours don't exist.
        % ask the user if he wants to run auto segmentation or if he knows
        % it will fail and we should just make some junky ones our of thin
        % air.
        
        selection = questdlg(['No Contours' get(handles.figure1,'Name') '?'],...
                             ['Contour files were not found.  Should I make them up? ' get(handles.figure1,'Name') '...'],...
                             'Yes','No','Yes');
        if strcmp(selection,'No')
            return;
        end
        theta = 0:.1:(2*pi);
        v.ContourData.AIF = horzcat((3*sin(theta)+v.sx/2)',(3*cos(theta)+v.sy/2)');
        v.ContourData.epi = horzcat((10*sin(theta)+v.sx/2)',(10*cos(theta)+v.sy/2)');
        v.ContourData.endo = horzcat((7*sin(theta)+v.sx/2)',(7*cos(theta)+v.sy/2)');
        v.ContourData.angle = pi;
        v.contours{1} = v.ContourData.AIF;
        v.contours{2} = v.ContourData.endo;
        v.contours{3} = v.ContourData.epi;
        for i=1:3
            v.masks(:,:,i) = roipoly(v.DCData.DCmri(:,:,1),v.contours{i}(:,1),v.contours{i}(:,2));
        end
    elseif(~UpdateOnlyRecent || temp.datenum > v.OldTimeStamp)
        v.ContourData.AIF = load(['Output/blood_polyCoords.study' num2str(v.seriesNum) '.slice' num2str(v.sliceNum)]);
        v.ContourData.epi = load(['Output/epi_polyCoords.study' num2str(v.seriesNum) '.slice' num2str(v.sliceNum)]);
        v.ContourData.endo = load(['Output/endo_polyCoords.study' num2str(v.seriesNum) '.slice' num2str(v.sliceNum)]);
        v.ContourData.angle = load(['Output/Roi_start_angle.study' num2str(v.seriesNum) '.slice' num2str(v.sliceNum)]);

        v.StartAngle = v.ContourData.angle;
        v.legacy = 0;
        if(v.ContourData.angle > 2*pi)
            v.legacy = 1;
            v.ContourData.angle = mod((v.ContourData.angle-90)*pi/180,2*pi);
        end
        v.contours = {};
        v.contours{1} = v.ContourData.AIF;
        v.contours{2} = v.ContourData.endo;
        v.contours{3} = v.ContourData.epi;
        v.center = ceil(mean(v.contours{3},1));
        for i=1:3
            v.masks(:,:,i) = roipoly(v.DCData.DCmri(:,:,1),v.contours{i}(:,1),v.contours{i}(:,2));
        end
    end
    
%     disp('In loadCR')
%     keyboard
    
    %set(handles.figure1,'UserData',v);
    sucessful = 1;

function KtransData = loadMostRecentKtrans(seriesNum,sliceNum)
    clear KtransData
    KtransData.ktrans = 0;
    files = dir(['Output/flowvalues.study' num2str(seriesNum) '.slice' num2str(sliceNum) '*']);
    times = zeros(length(files),1);
    for i=1:length(files)
        times(i) = files.datenum;
    end
    [times, IX] = sort(times); %#ok<ASGLU>
    files = files(IX);
    KtransData.ktrans = getKtrans(['Output/' files(1).name]);

% --- display stuff
function showImage(handles)
global v
    global showImageMutex
    if(~isfield(handles,'figure1')), return, end
    if(showImageMutex ~= 0), return; end
    showImageMutex = rand(1,1);
    history = showImageMutex;

    %v = get(handles.figure1,'UserData');
    %mycolormap = hsv;  %change this to get another colormap
    %colors = mycolormap(ceil(linspace(1,length(mycolormap),v.CurvesData.NAzimuthalRegions)),:);
    colors = lines(v.CurvesData.NAzimuthalRegions);

    %% main plot
    axes(handles.axes1); %#ok<MAXES>
    backgroundImage = getBackgroundImage(handles);
    if(v.showSmooth)
        pcolor(backgroundImage), shading interp, colormap gray
        if(~strcmp(v.BackgroundImageType,'Frames') && ~strcmp(v.BackgroundImageType,'Average'))
            colormap jet
        end
    else
        if(v.showMyoMask)
            mymin = min(min(backgroundImage));
            mymax = max(max(backgroundImage));
            temp = (backgroundImage -mymin)/(mymax-mymin);
            mask = (v.masks(:,:,3) - v.masks(:,:,2)) > 0;
            img(:,:,1) = temp;
            img(:,:,3) = temp;
            temp(mask) = sqrt(temp(mask));
            img(:,:,2) = temp;
            imagesc(img), colormap gray
        else
            imagesc(backgroundImage), colormap gray
            if(~strcmp(v.BackgroundImageType,'Frames') && ~strcmp(v.BackgroundImageType,'Average'))
                colormap jet
            end
        end
    end
    if(history ~= showImageMutex), 
        %v = get(handles.figure1,'UserData'); 
    end
    axis ij
    hold on

    %contours
    %make a smooth version of the curves
    blood = v.contours{1};
    endo = v.contours{2};
    epi = v.contours{3};

    v.showSmooth=0; 
    if(v.showSmooth)
        endo = vertcat(endo((end-10):end,:),endo,endo(1:10,:));
        epi = vertcat(epi((end-10):end,:),epi,epi(1:10,:));
        blood = vertcat(blood((end-5):end,:),blood,blood(1:5,:));
        endo(:,1) = filter2(fspecial('gaussian',size(endo(:,1)),2),endo(:,1));
        endo(:,2) = filter2(fspecial('gaussian',size(endo(:,2)),2),endo(:,2));
        epi(:,1) = filter2(fspecial('gaussian',size(epi(:,1)),2),epi(:,1));
        epi(:,2) = filter2(fspecial('gaussian',size(epi(:,2)),2),epi(:,2));
        blood(:,1) = filter2(fspecial('gaussian',size(blood(:,1)),1),blood(:,1));
        blood(:,2) = filter2(fspecial('gaussian',size(blood(:,2)),1),blood(:,2));
        endo = endo(11:(end-11),:);
        epi = epi(11:(end-11),:);
        blood = blood(6:(end-6),:);
    end
    type = '-';
    if(v.DeletePressed && v.contourSelected == 1)
        type = [type '-'];
    end
    plot([blood(:,1)' blood(1,1)],[blood(:,2)' blood(1,2)],'k','LineWidth',(v.contourSelected == 1)+1,'LineStyle',type);
    center = mean(endo,1);

    %if(history ~= showImageMutex), v = get(handles.figure1,'UserData'); end

    [endotheta,endoR] = cart2pol(endo(:,2) - center(2),endo(:,1) - center(1));
    [epitheta,~] = cart2pol(epi(:,2) - center(2),epi(:,1) - center(1));
    
    regionalEpi = mod(epitheta - v.StartAngle,2*pi);
    labeledepi = floor(regionalEpi/(2*pi)*v.CurvesData.NAzimuthalRegions);
    if(~v.legacy)
        labeledepi = mod(labeledepi+2,6);
    end
    
    regionalEndo = mod(endotheta - v.StartAngle,2*pi);
    labeledendo = floor(regionalEndo/(2*pi)*v.CurvesData.NAzimuthalRegions);
    if(~v.legacy)
        labeledendo = mod(labeledendo+2,6);
    end

    for i=1:v.CurvesData.NAzimuthalRegions
        mask = (labeledendo+1)==i;
        endot = squeeze(endo(mask,:));
        % extend mask by 1 in both directions
       % find(diff(mask))=0)
       ind1=find(diff(mask)==1);
       %mask(ind1-1)=1;
       ind2=find(diff(mask)==-1);
       %mask(ind2)=1;
        
        endot = squeeze(endo(mask,:));
     %   keyboard
        timesShifted=0;
       try  while(max(abs(diff(endot(:,1)))) > 3 || max(abs(diff(endot(:,2)))) > 3)
            if(timesShifted > length(endot))
                break;
            end
            timesShifted = timesShifted + 1;
            endot = circshift(endot,[1 0]);
           end
       catch
           disp(' draw region 2, endo ')
           selection = questdlg(['Problem with Contours likely from 3.1, seems not enough points.  Want some quick junky ones that can then be redrawn with 3.11? ' get(handles.figure1,'Name') '?'],...
                             ['Contour files ' get(handles.figure1,'Name') '...'],...
                             'Yes','No','Yes');
        if strcmp(selection,'No')
            return;
        end
        theta = 0:.1:(2*pi);
        v.ContourData.AIF = horzcat((3*sin(theta)+v.sx/2)',(3*cos(theta)+v.sy/2)');
        v.ContourData.epi = horzcat((10*sin(theta)+v.sx/2)',(10*cos(theta)+v.sy/2)');
        v.ContourData.endo = horzcat((7*sin(theta)+v.sx/2)',(7*cos(theta)+v.sy/2)');
        v.ContourData.angle = pi;
        v.contours{1} = v.ContourData.AIF;
        v.contours{2} = v.ContourData.endo;
        v.contours{3} = v.ContourData.epi;
        for i=1:3
            v.masks(:,:,i) = roipoly(v.DCData.DCmri(:,:,1),v.contours{i}(:,1),v.contours{i}(:,2));
        end
        
        
        for i=1:3
        v.contourSelected = i;
        %set(handles.figure1,'UserData',v);
        MasksChanged(handles);
        %v = get(handles.figure1,'UserData');
    end

%     %set(handles.figure1,'UserData',v);
%     set(handles.radiobutton3,'Value',v.contourSelected == 3);
%     set(handles.radiobutton4,'Value',v.contourSelected == 2);
%     set(handles.radiobutton5,'Value',v.contourSelected == 1);
%     
%     showImage(handles);
%     showCurves(handles);
%     showShifts(handles);
%     
    
    
    break;
        
        
        
%            v.DeletePressed=1   % if not enough points drawn (like from 3.1, force redraw here
%            v.contourSelected = 2
%            showImage(handles)
%         v.DeletePressed = 0;
%         [v.masks(:,:,v.contourSelected),X,Y] = roipoly();
%         if(v.contourSelected==2)  % endo? then dilate
%             se=strel('disk',1);
%             img_dilated=imdilate(v.masks(:,:,v.contourSelected), se);
%             [r,c] = find(img_dilated); 
%             v.contours{v.contourSelected}=bwtraceboundary(img_dilated,[r(1) c(1)],'N');
%             MasksChanged(handles,1:v.st,0);
%         else
%             v.contours{v.contourSelected} = horzcat(X,Y);
%         %set(handles.figure1,'UserData',v);
%             MasksChanged(handles,1:v.st,1);  %EVRD, odd, freezes otherewise if endo or epi contours.. 
%         %  dilate/erode won't work after do a delete,
%         end
%         showImage(handles)
%         set(handles.figure1,'WindowScrollWheelFcn',@(hObject,eventdata)CR_editor('figure1_WindowScrollWheelFcn',hObject,eventdata,guidata(hObject)));
%     
        
        
       end
       
        type = '-';
        if(v.DeletePressed && v.contourSelected == 2)
            type = [type '-']; %#ok<AGROW>
        end
        plot(endot(:,1), endot(:,2), 'Color',colors(i,:),'LineWidth',(v.contourSelected == 2)+1,'LineStyle',type);
%         x1poly=endot(:,1); y1poly=endot(:,2);
% h1=line(x1poly(1:end-1),y1poly(1:end-1));
% set(h1,'Color',colors(i,:));  %[0 1 0]); 
% set(h1,'Linewidth',(v.contourSelected == 2)+1)
% set(h1,'LineStyle',type)
        
        mask = (labeledepi+1) == i;
        epit = squeeze(epi(mask,:));
        if(length(epit) < 2)
            continue;
        end
        timesShifted = 0;
       try
           while(max(abs(diff(epit(:,1)))) > 3 || max(abs(diff(epit(:,2)))) > 3)
            if(timesShifted > length(epit))
                break;  %return;
            end
            timesShifted = timesShifted + 1;
            epit = circshift(epit,[1 0]);
        end
        
        catch
           disp('choose region 3, epi ')
       end
    
        type = '-';
        if(v.DeletePressed && v.contourSelected == 3)
            type = [type '-']; %#ok<AGROW>
        end
        plot(epit(:,1), epit(:,2), 'Color',colors(i,:),'LineWidth',(v.contourSelected == 3)+1,'LineStyle',type);
       % imagesc(backgroundImage), colormap jet

        %if(history ~= showImageMutex), v = get(handles.figure1,'UserData'); end

        %k-trans values overlays
        temp = epitheta(mask);
        if(any(temp < -3*pi/4) && any(temp > 3*pi/4))
            temp(temp<0) = temp(temp<0) + 2*pi;
        end
        textAngle = median(temp);
        textRadius = 1.8*max(endoR);
        [textPositionY, textPositionX] = pol2cart(textAngle,textRadius);
        textPositionX = textPositionX+center(1);
        textPositionY = textPositionY+center(2);
        text(textPositionX,textPositionY,{...
            ['Reg. ' num2str(i) ' : ' num2str(v.KtransData.ktrans(i))],...
            },'Color',colors(i,:));
    end
%          x1poly=epi(:,1); y1poly=epi(:,2);
% %bw1=roipoly(tmpImg,x1poly,y1poly);
% h1=line(x1poly,y1poly);
% set(h1,'Color',[0 1 0]); 
% set(h1,'Linewidth',0.5)
% x1poly=endo(:,1); y1poly=endo(:,2);
% h1=line(x1poly,y1poly);
% set(h1,'Color',[0 1 1]); 
% set(h1,'Linewidth',0.5)
%curves = ExtractCurves(cinemri1,blood,endo,epi,angle,numAzimuthalRegions,referenceFrame)



    %show the starting angle
    [AnglePointy,AnglePointx] = pol2cart(repmat(v.StartAngle,1,length(10:99)),10:99);
    AnglePointx = AnglePointx + center(1);AnglePointy = AnglePointy + center(2);
    if(sum(AnglePointx>10) > sum(AnglePointy>10))
        AnglePointx = AnglePointx(AnglePointx>10);AnglePointy = AnglePointy(AnglePointx>10);
    else
        AnglePointx = AnglePointx(AnglePointy>10);AnglePointy = AnglePointy(AnglePointy>10);
    end
    plot(AnglePointx,AnglePointy);


     %v = get(handles.figure1,'UserData'); 

    % display circle for bumping tool
    if(v.bumping && v.mouse(1) ~= 1 && v.mouse(2) ~= 1 && v.mouse(2) < size(v.masks,1) && v.mouse(1) < size(v.masks,2))
        bumpColors = 'br';
        symbols = {'-','+'};
        inOut = v.masks(ceil(v.mouse(2)),ceil(v.mouse(1)),v.contourSelected);
        plot(v.r*sin(0:.1:(2*pi)) + v.mouse(1),v.r*cos(0:.1:(2*pi)) + v.mouse(2),'Color',bumpColors(inOut+1));
        text(v.mouse(1),v.mouse(2),symbols{inOut+1},'Color',bumpColors((~inOut)+1));
    end
    
    showImageMutex = 0;
    hold off
function img = getBackgroundImage(handles)
global v
if(strcmp(v.BackgroundImageType,'Frames'))
    img = v.DCData.DCmri(:,:,v.t);
elseif(strcmp(v.BackgroundImageType,'Variance'))
    img = var(v.DCData.DCmri,0,3);
elseif(strcmp(v.BackgroundImageType,'Standard Deviation'))
    img = std(v.DCData.DCmri,0,3);
elseif(strcmp(v.BackgroundImageType,'Average'))
    img = mean(v.DCData.DCmri,3);
elseif(strcmp(v.BackgroundImageType,'Delays') || strcmp(v.BackgroundImageType,'Upslope'))
    if(~isfield(v,'Upslope'))
        myslopefinder = -diff(fspecial('gaussian',[(v.st+1),1],7));
        
        AIF_upslope = max(filter2(myslopefinder,v.blood));
        
        v.Upslope = zeros(v.sx,v.sy);
        v.Delays = zeros(v.sx,v.sy);
        for x=1:v.sx
            for y=1:v.sy
                [v.Upslope(x,y) v.Delays(x,y)] = max(filter2(myslopefinder,squeeze(v.DCData.DCmri(x,y,:))));
            end
        end
        v.Upslope = v.Upslope ./ AIF_upslope;
        mystd = std(v.Delays(:));
        mymean = mean(v.Delays(:));
        v.Delays(v.Delays-mymean> 2*mystd) = 0;
        v.Delays(v.Delays-mymean< -2*mystd) = 0;
        %set(handles.figure1,'UserData',v);
    end
    if(strcmp(v.BackgroundImageType,'Delays'))
        img = v.Delays;
    else
        img = v.Upslope;
    end
elseif(strcmp(v.BackgroundImageType,'Manhattan Distance'))
    if(~isfield(v,'Manhattan'))
        tissuCurve = mean(v.CurvesData.curves,1)';
        v.Manhattan = zeros(v.sx,v.sy);
        for x=1:v.sx
            for y=1:v.sy
                v.Manhattan(x,y) = norm(squeeze(v.DCData.DCmri(x,y,:)) - tissuCurve);
            end
        end
        %set(handles.figure1,'UserData',v);
    end
    img = v.Manhattan;
elseif(...
        strcmp(v.BackgroundImageType,'euclidean') || ...
        strcmp(v.BackgroundImageType,'seuclidean') || ...
        strcmp(v.BackgroundImageType,'cityblock') || ...
        strcmp(v.BackgroundImageType,'minkowski') || ...
        strcmp(v.BackgroundImageType,'chebychev') || ...
        strcmp(v.BackgroundImageType,'mahalanobis') || ...
        strcmp(v.BackgroundImageType,'cosine') || ...
        strcmp(v.BackgroundImageType,'correlation') || ...
        strcmp(v.BackgroundImageType,'spearman')...
        )
    if(~isfield(v,v.BackgroundImageType))
        tissuCurve = mean(v.CurvesData.curves,1);
        for x=1:v.sx
            v.(v.BackgroundImageType)(x,:) = pdist2(squeeze(v.DCData.DCmri(x,:,:)),tissuCurve,v.BackgroundImageType);
        end
    end
    img = v.(v.BackgroundImageType);
else
    img = v.DCData.DCmri(:,:,v.t);
end

function showCurves(handles)
global v
    %% curves

    if(~isfield(handles,'figure1')), return, end
    %v = get(handles.figure1,'UserData');
    colors = lines(v.CurvesData.NAzimuthalRegions);
    % circ shift by two
    %colors = circshift(colors,2);
    axes(handles.axes2);cla; %#ok<MAXES>
    hold on
    if(v.contourSelected == 1)
        plot(v.blood,'k');
    end
    for i=1:v.CurvesData.NAzimuthalRegions
        plot(v.CurvesData.deltaSIcurves(i,:),'Color',colors(i,:));
    end
    plot([v.t v.t],[min(v.CurvesData.deltaSIcurves(:)) max(max(v.CurvesData.deltaSIcurves(2:end,:)))],'k'); hold off
    hold off

function showShifts(handles)
global v
    %% show the shifts
    if(~isfield(handles,'figure1')), return, end
    %v = get(handles.figure1,'UserData');
    if(~isfield(v,'ShiftsData'))
        return;
    end
    myshifts = zeros(length(v.ShiftsData),2);
    for t=1:length(v.ShiftsData)
        myshifts(t,:) = v.ShiftsData(t).shift;
    end

    axes(handles.axes3); %#ok<MAXES>
    
    
    plot(myshifts(:,1));
    hold on, plot([v.t v.t],[min(min(myshifts(:,1)),-1) max(max(myshifts(:,1)),1)],'k'); hold off

    set(gca,'XTick',[]);
    set(gca,'YAxisLocation','right');
    axes(handles.axes4) %#ok<MAXES>
    
    plot(myshifts(:,2));

    hold on, plot([v.t v.t],[min(min(myshifts(:,2)),-1) max(max(myshifts(:,2)),1)],'k'); hold off
    set(gca,'XTick',[]);
    set(gca,'YAxisLocation','right');



function MasksChanged(handles,t,dilate_or_erode)
global v MasksChangedMutex
    %% --- updates the contours
    if(~isfield(handles,'figure1')), return, end

    if(exist('t','var'))
        ranget = t;
    else
        ranget = 1:v.st;
    end
    
    MasksChangedMutex = rand(1,1);
    myhistory = MasksChangedMutex;
    %v = get(handles.figure1,'UserData');
    myo = (v.masks(:,:,3) - v.masks(:,:,2)) > 0;
    [row,col] = find(myo);
    v.center = ceil(mean(v.contours{3},1));
    [theta, ~] = cart2pol(row - v.center(2), col - v.center(1));

    blood = v.masks(:,:,1);
    for t=ranget
        temp = v.DCData.DCmri(:,:,t);
        v.blood(t) = mean(temp(blood>0));
    end
    v.blood = v.blood - mean(v.blood(1:5));
    
    theta = mod(theta - v.StartAngle,2*pi);
    labeled = floor(theta/(2*pi)*v.CurvesData.NAzimuthalRegions);
    if(~v.legacy)
        labeled = mod(labeled+2,6);
    end

    referenceFrame=25;
    curves = ExtractCurves_useWith_CReditor(v.DCData.DCmri,v.contours{1},v.contours{2},v.contours{3},v.StartAngle,v.CurvesData.NAzimuthalRegions,referenceFrame);
    v.blood= curves(1,ranget);
    v.blood= v.blood - mean(v.blood(1:5));
    curves=curves(2:end,:);  %skip AIF, gets shown only if AIF contour selected
    % now convert to deltaSI
    for i = 1:v.CurvesData.NAzimuthalRegions
          v.CurvesData.deltaSIcurves(i,ranget) = curves(i,ranget) - median(curves(i,1:5));
    end

    
%     for i = 1:v.CurvesData.NAzimuthalRegions
%         %find what region each point in the myo is in
%         maskedrow = row(labeled+1 == i);
%         maskedcol = col(labeled+1 == i);
%         if(isempty(maskedrow))
%             disp(['Empty region' num2str(i+1)]);
%             continue;
%         end
%         %effectivly do pixelwise
%         clear temp
%         for j=1:length(maskedrow)
%             temp(j,ranget) = v.DCData.DCmri(maskedrow(j),maskedcol(j),ranget);
%         end
%         if(MasksChangedMutex ~= myhistory)
%             v.CurvesData.deltaSIcurves(i,ranget) = mean(temp(randi(size(temp,1),min(size(temp,2),100)),ranget),1);
%         else
%             v.CurvesData.deltaSIcurves(i,ranget) = mean(temp,1);
%         end
%         v.CurvesData.deltaSIcurves(i,ranget) = v.CurvesData.deltaSIcurves(i,ranget) - median(v.CurvesData.deltaSIcurves(i,1:5));
%     end
% 
%     
    
    
    
    
    %fix the contours
    [r,c] = find(v.masks(:,:,v.contourSelected));  %maybe the AIF does not overlap the endocardium center
    contour = bwtraceboundary(v.masks(:,:,v.contourSelected),[r(1) c(1)],'N');    
    %EVRD 12/20/10
    if(exist('dilate_or_erode','var'))   % only time will go from mask to create contours. Otherwise get burned as shrinks!
     v.contours{v.contourSelected} = contour(:,[2 1]);
     hold on; plot(contour(:,2),contour(:,1),'g','LineWidth',1)
    end   
    
    %from Brian 12/20/10  -didn't just work when tested, EVRD, shrunk AIF, maybe offcenter too.  maybe
    %revisit someday
%     %do polar rounding
%     [Theta,R] = cart2pol(contour(:,1)-v.center(1),contour(:,2)-v.center(2));
%     R = R*1.01;
%     [contour(:,1),contour(:,2)] = pol2cart(Theta,R);
%     contour(:,1) = contour(:,1)+v.center(1);
%     contour(:,2) = contour(:,2)+v.center(2);
%     contour = contour(:,[2 1]);
%     %save it away
     %v.contours{v.contourSelected} = contour(:,[2 1]);;
% hold on; plot(contour(:,2),contour(:,1),'g','LineWidth',1)

    %set(handles.figure1,'UserData',v);
    showCurves(handles);

% --- Outputs from this function are returned to the command line.
function varargout = CR_editor_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                         ['Close ' get(handles.figure1,'Name') '...'],...
                         'Yes','No','Yes');
    if strcmp(selection,'No')
        return;
    end
    saveRC(handles,1);
    delete(handles.figure1)

function saveRC(handles,ask)
global v
    %v = get(handles.figure1,'UserData');
    if(~exist('ask','var'))
        ask = 1;
    end
    if(~isfield(v,'DCData'))
        return;
    end
    
    if(ask)
        selection = questdlg('Save Movie?',...
                            'Save things',...
                            'Yes to All','Yes','No to All','Yes');
    end
    if((ask && strcmp(selection,'Yes to All')))
        ask = 0;
    end
    if((ask && strcmp(selection,'No to All')))
        return;
    end
    if(~ask || (ask && strcmp(selection,'Yes')))
        disp('Saving Movie');
        cinemri1 = v.DCData.DCmri; %#ok<NASGU>
        save(['Output/cinemri1.study' num2str(v.seriesNum) '.slice' num2str(v.sliceNum) '.mat'],'cinemri1');
    else
        disp('Not movie');
    end
    if(v.changedShifts)
        if(ask)
            selection = questdlg('Save Shifts?',...
                                'Shifts',...
                                'Yes','No','Yes');
        end
        if(~ask || (ask && strcmp(selection,'Yes')))
            disp('Saving Manual shifts');
            shifts = zeros(v.st,2);
            for t=1:v.st
                shifts(t,:) = v.ShiftsData(t).shift;
            end
            save(['Output/shiftsMAN.study' num2str(v.seriesNum) '.slice' num2str(v.sliceNum) '.txt'],'shifts','-ascii','-tabs');
        else
            disp('Not saving shifts');
        end
    end
    if(1 || v.changedContours ~= 0)
        if(ask)
            selection = questdlg('Save Contours?',...
                                'Contours',...
                                'Yes','No','Yes');
        end
        if(~ask || (ask && strcmp(selection,'Yes')))
            epi = v.contours{3}; %#ok<NASGU>
            endo = v.contours{2}; %#ok<NASGU>
            AIF = v.contours{1}; %#ok<NASGU>
            angle = v.StartAngle; %#ok<NASGU>
            save(['Output/blood_polyCoords.study' num2str(v.seriesNum) '.slice' num2str(v.sliceNum)],'AIF','-ascii','-double');
            save(['Output/epi_polyCoords.study' num2str(v.seriesNum) '.slice' num2str(v.sliceNum)],'epi','-ascii','-double');
            save(['Output/endo_polyCoords.study' num2str(v.seriesNum) '.slice' num2str(v.sliceNum)],'endo','-ascii','-double');
            save(['Output/Roi_start_angle.study' num2str(v.seriesNum) '.slice' num2str(v.sliceNum)],'angle','-ascii','-double');    
            curvefilename=strcat('Output/curves.study',num2str(v.seriesNum),'.slice',num2str(v.sliceNum),'.mat');   %EVRD 12/28/10, not sure where to put this, want to
            curves=[v.blood; v.CurvesData.curves];  % write out if contours or shifts are changed... 
            save(curvefilename, 'curves')
        else
            disp('Not saving contours');
        end
    end
%     disp('In saveRC')
%     keyboard
return


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
global v
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%v = get(handles.figure1,'UserData');
v.contourSelected = 3;
%set(handles.figure1,'UserData',v);
showImage(handles)

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
global v
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%v = get(handles.figure1,'UserData');
v.contourSelected = 2;
%set(handles.figure1,'UserData',v);
showImage(handles)

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
global v
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%v = get(handles.figure1,'UserData');
v.contourSelected = 1;
%set(handles.figure1,'UserData',v);
showImage(handles)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%v = get(handles.figure1,'UserData');
%objects = findobj(gcf,'style','push');
global v
if(v.bumping)
    v.bumping = 0;
    set(handles.pushbutton6,'FontWeight','normal');
%     %restore the callbacks of all controls
%     for i=1:length(objects)
%         set(objects(i), 'keypressfcn',v.keyPressFunctionBackup(objects(i)));
%     end
else
    v.bumping = 1;
    set(handles.pushbutton6,'FontWeight','Bold');
%     %save away the keydown callbacks of all controls
%     if(~isfield(v,'keyPressFunctionBackup'))
%         v.keyPressFunctionBackup = containers.Map({objects(1)},{{get(objects(1),'keypressfcn')}});, remove(v.keyPressFunctionBackup,objects(1));
%     end
%     for i=1:length(objects)
%         v.keyPressFunctionBackup(objects(i)) = get(objects(i),'keypressfcn');
%     end
%     set(objects,'keypressfcn',@figure1_WindowKeyPressFcn)
end
%set(handles.figure1,'UserData',v);
showImage(handles);





% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
global v
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%v = get(handles.figure1,'UserData');
if(v.bumping)
    pushbutton6_Callback(hObject, eventdata, handles)
end
%v = get(handles.figure1,'UserData');
v.waitingForAngleClick = 1;
%set(handles.figure1,'UserData',v);

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
global v
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
global keyPressMutex
    if(keyPressMutex ~= 0), return; end
    keyPressMutex = 1; %#ok<NASGU>
    if(~exist('handles','var'))
        keyPressMutex = 0;
        return;
    end
    %disp('------------------');
    objects = fieldnames(handles);

    for i=1:length(objects)
        if(hObject == handles.(objects{i}))
            %disp(['I got called by ' objects{i}]);
        end
    end
    %v = get(handles.figure1,'UserData');
    v.HitAKeyTimeStamp = now;
    %set(handles.figure1,'UserData',v);
    if(~isfield(eventdata,'Key'))
        keyPressMutex = 0;
        return;   % for some reason this function is getting called from scroll
    end
    if(strcmp(eventdata.Key,'p') || strcmp(eventdata.Key,'space'))
        if(v.playing)
            v.playing = 0;
            %set(handles.figure1,'UserData',v);
        else
            v.playing = 1;
            %set(handles.figure1,'UserData',v);
            keyPressMutex = 0; %#ok<NASGU>
            while(1)
                %v = get(handles.figure1,'UserData');
                if(v.playing==0), break; end
                v.t = v.t+1;
                if(v.t > v.st)
                    v.t = 1;
                end
                %set(handles.figure1,'UserData',v);
                showImage(handles);
                showShifts(handles)
                showCurves(handles)
                drawnow
            end
        end
    elseif(strcmp(eventdata.Key,'1'))
        v.contourSelected = 1;
        set(handles.radiobutton3,'Value',v.contourSelected == 3);
        set(handles.radiobutton4,'Value',v.contourSelected == 2);
        set(handles.radiobutton5,'Value',v.contourSelected == 1);
        %set(handles.figure1,'UserData',v);
        showImage(handles)
    elseif(strcmp(eventdata.Key,'2'))
        v.contourSelected = 2;
        set(handles.radiobutton3,'Value',v.contourSelected == 3);
        set(handles.radiobutton4,'Value',v.contourSelected == 2);
        set(handles.radiobutton5,'Value',v.contourSelected == 1);
        %set(handles.figure1,'UserData',v);
        showImage(handles)
    elseif(strcmp(eventdata.Key,'3'))
        v.contourSelected = 3;
        set(handles.radiobutton3,'Value',v.contourSelected == 3);
        set(handles.radiobutton4,'Value',v.contourSelected == 2);
        set(handles.radiobutton5,'Value',v.contourSelected == 1);
        %set(handles.figure1,'UserData',v);
        showImage(handles)
    elseif(strcmp(eventdata.Key,'uparrow') || ...
            strcmp(eventdata.Key,'downarrow') || ...
            strcmp(eventdata.Key,'leftarrow') || ...
            strcmp(eventdata.Key,'rightarrow'))
        if(strcmp(eventdata.Key,'uparrow')), shift = [-1 0]; end
        if(strcmp(eventdata.Key,'downarrow')), shift = [1 0]; end
        if(strcmp(eventdata.Key,'leftarrow')), shift = [0 -1]; end
        if(strcmp(eventdata.Key,'rightarrow')), shift = [0 1]; end
       % v.DCData.DCmri = circshift(v.DCData.DCmri,shift);
        %get the timestamp the user is looking at
        timestamp = v.DCData.timestamp(v.t);
        v.changedShifts = 1;
        temp = [v.ShiftsData(:).timestamp];
        index = find(temp == timestamp);
        v.ShiftsData(index).shift = v.ShiftsData(index).shift + shift;
        v.DCData.DCmri(:,:,index) = circshift(v.DCData.DCmri(:,:,index),shift);  %EVRD  12/14/10
       % set(handles.figure1,'UserData',v);
        MasksChanged(handles)
        showImage(handles)
        showShifts(handles)
    elseif(strcmp(eventdata.Key,'delete'))
        v.DeletePressed = 1;
        %set(handles.figure1,'UserData',v);
        showImage(handles)
        v.DeletePressed = 0;
        [v.masks(:,:,v.contourSelected),X,Y] = roipoly();
        if(v.contourSelected==2)  % endo? then dilate
            se=strel('disk',1);
            img_dilated=imdilate(v.masks(:,:,v.contourSelected), se);
            [r,c] = find(img_dilated); 
            v.contours{v.contourSelected}=bwtraceboundary(img_dilated,[r(1) c(1)],'N');
            MasksChanged(handles,1:v.st,0);
        else
            v.contours{v.contourSelected} = horzcat(X,Y);
        %set(handles.figure1,'UserData',v);
            MasksChanged(handles,1:v.st,1);  %EVRD, odd, freezes otherewise if endo or epi contours.. 
        %  dilate/erode won't work after do a delete,
        end
        showImage(handles)
        set(handles.figure1,'WindowScrollWheelFcn',@(hObject,eventdata)CR_editor('figure1_WindowScrollWheelFcn',hObject,eventdata,guidata(hObject)));
    elseif(strcmp(eventdata.Key,'pagedown'))
        v.t = v.t+1;
        if(v.t > size(v.DCData.DCmri,3))
            v.t = 1;
        end
        %set(handles.figure1,'UserData',v);
        showImage(handles)
        showCurves(handles)
        showShifts(handles)
    elseif(strcmp(eventdata.Key,'pageup'))
        v.t = v.t-1;
        if(v.t ==0)
            v.t = size(v.DCData.DCmri,3);
        end
        %set(handles.figure1,'UserData',v);
        showImage(handles)
        showCurves(handles)
        showShifts(handles)
    elseif(strcmp(eventdata.Key,'escape'))
        if(v.bumping)
            %set(handles.figure1,'UserData',v);
            pushbutton6_Callback(hObject, eventdata, handles);
            showImage(handles)
        else
            CloseMenuItem_Callback(hObject, eventdata, handles);
        end
        keyPressMutex = 0;
        return;
    elseif(strcmp(eventdata.Key,'equal') || strcmp(eventdata.Key,'plus'))
        erodeDilate(handles,1);
        showImage(handles)
    elseif(strcmp(eventdata.Key,'hyphen') || strcmp(eventdata.Key,'minus'))
        erodeDilate(handles,-1);
        showImage(handles)
    end
    
    keyPressMutex = 0;

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global v
    %v = get(handles.figure1,'UserData');
    v.showSmooth = ~v.showSmooth;
    onoff = {'off','on'};
    set(hObject,'Checked',onoff{v.showSmooth+1});
    %set(handles.figure1,'UserData',v);
    showImage(handles);


% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dragMutex v
    if(dragMutex ~= 0)
        return; 
    end
    dragMutex = 1; %#ok<NASGU>
    %v = get(handles.figure1,'UserData');
    mouse = get(handles.axes1,'CurrentPoint');
    mouse = mouse(1,1:2);
    changedSomething = 0;
    %find what axes the user is hovering over
    if(~isfield(v,'masks')), 
        dragMutex = 0;
        return; 
    end;
    if(mouse(2) < size(v.masks,1) && mouse(1) < size(v.masks,2))
        if(v.bumping)
            changedSomething = 1;
        end
        changedMask = 0;
        v.mouse = mouse;
        if(v.waitingForAngleClick)
            [theta, ~] = cart2pol(v.mouse(1) - v.center(1), v.mouse(2) - v.center(2));
            v.StartAngle = -theta+pi/2;
            changedSomething = 1;
        end
        if(strcmp(v.clickType,'normal') && v.bumping)
            mask = v.masks(:,:,v.contourSelected);
            inOut = mask(ceil(v.mouse(2)),ceil(v.mouse(1)));
            [X, Y] = meshgrid(1:(2*v.r+1), 1:(2*v.r+1));
            rangex = max(1,ceil(v.mouse(2))-v.r) :min(size(v.masks,1),ceil(v.mouse(2))+v.r);
            rangey = max(1,ceil(v.mouse(1))-v.r) :min(size(v.masks,2),ceil(v.mouse(1))+v.r);
            indexMask = (sqrt((X-v.r-1).^2 + (Y-v.r-1).^2) <= v.r);
            temp = mask(rangex, rangey);
            if(inOut)
                temp = ~temp;
            end
            if(all(size(temp) == size(indexMask)) && any(temp(indexMask)) > 0)
                changedMask = 1;
                temp(indexMask) = 0;
                if(inOut)
                     temp = ~temp;
                end
                mask(rangex, rangey) = temp;
                v.masks(:,:,v.contourSelected) = mask;
            end

        end
        if(strcmp(v.clickType,'normal') && ~v.bumping)
            if((v.mouse(1) == 1 && v.mouse(2) == 1) || (v.oldmouse(1) == 1 && v.oldmouse(2) == 1))
                v.oldmouse = v.mouse;
                %set(handles.figure1,'UserData',v);
                showImage(handles);
                dragMutex = 0;
                return;
            end
            if(norm(v.oldmouse - v.mouse) > 1)
                diff = round(v.oldmouse - v.mouse);
                diff = -diff([2 1]);
                v.changedShifts = 1;
                changedSomething = 1;
                v.oldmouse = v.mouse;
                v.ShiftsData(v.t).shift = v.ShiftsData(v.t).shift + diff;
                v.DCData.DCmri(:,:,v.t) = circshift(v.DCData.DCmri(:,:,v.t),diff);
                %set(handles.figure1,'UserData',v);
                
                showCurves(handles);
                showShifts(handles);
            end
        end
        if(strcmp(v.clickType,'alt'))
            if((v.mouse(1) == 1 && v.mouse(2) == 1) || (v.oldmouse(1) == 1 && v.oldmouse(2) == 1))
                v.oldmouse = v.mouse;
                %set(handles.figure1,'UserData',v);
                dragMutex = 0;
                return;
            end
            if(norm(v.oldmouse - v.mouse) > 1)
                diff = round(v.oldmouse - v.mouse);
                v.oldmouse = v.mouse;
                v.masks(:,:,v.contourSelected) = circshift(v.masks(:,:,v.contourSelected),-diff([2 1]));
                changedMask = 1;
                changedSomething = 1;
            end
        end
    else
        v.mouse = [1 1];
    end
    if(norm(v.oldmouse - v.mouse) > 1)
        v.oldmouse = v.mouse;
    end

    if(changedSomething)
        %set(handles.figure1,'UserData',v);
        if(changedMask)
            v.changedContours = 1;
            MasksChanged(handles);
        end
        showImage(handles);
    end
    dragMutex = 0;

% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
global v
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)
    if(~isfield(handles,'figure1')), return, end
    %v = get(handles.figure1,'UserData');
    if(v.bumping)
        v.r = v.r+eventdata.VerticalScrollCount;
        if(v.r < 0), v.r = 1;end
        %set(handles.figure1,'UserData',v);
    else
        erodeDilate(handles,eventdata.VerticalScrollCount);
    end
    showImage(handles);

function erodeDilate(handles,count)
    global v
        directions = {'erode','dilate'};
        v.masks(:,:,v.contourSelected) = bwmorph(v.masks(:,:,v.contourSelected),directions{(count>0)+1},abs(count));
        if(sum(sum(v.masks(:,:,v.contourSelected))) <8)
            v.masks((v.center(2)-1):(v.center(2)+1),(v.center(1)-1):(v.center(1)+1),v.contourSelected) = 1;
        end
        v.changedContours = 1;
        %set(handles.figure1,'UserData',v);
        dilate_or_erode=1;
        t=1:v.st;  % just passing in default, so will see dilated
        MasksChanged(handles,t,dilate_or_erode);
% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
global v
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    if(~isfield(handles,'figure1')), return, end
    %v = get(handles.figure1,'UserData');
    v.clickType = get(handles.figure1,'SelectionType');
    if(v.waitingForAngleClick)
        v.waitingForAngleClick = 0;
        [theta, ~] = cart2pol(v.mouse(1) - v.center(1), v.mouse(2) - v.center(2));
        v.StartAngle = -theta+pi/2;
        v.changedContours = 1;
        %set(handles.figure1,'UserData',v);
        MasksChanged(handles);
    else
        %set(handles.figure1,'UserData',v);
    end

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
global v
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    %v = get(handles.figure1,'UserData');
    %set(handles.figure1,'UserData',v);
if(strcmp(v.clickType,'normal') && ~v.bumping)
    
    MasksChanged(handles);
end
    v.clickType = '';


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu7.
function popupmenu7_Callback(hObject, eventdata, handles)
global v
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu7 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu7
    %v = get(handles.figure1,'UserData');
    contents = cellstr(get(hObject,'String'));
    v.BackgroundImageType = contents{get(hObject,'Value')};
    %set(handles.figure1,'UserData',v);
    showImage(handles)

% --- Executes during object creation, after setting all properties.
function popupmenu7_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles) %#ok<*INUSL,*DEFNU>
global v
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
    %v = get(handles.figure1,'UserData');
    v.contourSelected = 3;
    set(handles.radiobutton3,'Value',1);
    set(handles.radiobutton4,'Value',0);
    set(handles.radiobutton5,'Value',0);
    %set(handles.figure1,'UserData',v);
    showImage(handles)
    showCurves(handles)

% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
global v
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4
    %v = get(handles.figure1,'UserData');
    v.contourSelected = 2;
    set(handles.radiobutton3,'Value',0);
    set(handles.radiobutton4,'Value',1);
    set(handles.radiobutton5,'Value',0);
    %set(handles.figure1,'UserData',v);
    showImage(handles)
    showCurves(handles)

% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
global v
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5
    %v = get(handles.figure1,'UserData');
    v.contourSelected = 1;
    set(handles.radiobutton3,'Value',0);
    set(handles.radiobutton4,'Value',0);
    set(handles.radiobutton5,'Value',1);
    %set(handles.figure1,'UserData',v);
    showImage(handles)
    showCurves(handles)


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
global v
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    %v = get(handles.figure1,'UserData');
    v.showMyoMask = ~v.showMyoMask;
    if(v.showMyoMask && v.showSmooth)
        %set(handles.figure1,'UserData',v);
        Untitled_2_Callback(hObject, eventdata, handles);
        %v = get(handles.figure1,'UserData');
    end
    onoff = {'off','on'};
    set(hObject,'Checked',onoff{v.showMyoMask+1});
    %set(handles.figure1,'UserData',v);
    showImage(handles);


% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
global v
% hObject    handle to Untitled_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    %v = get(handles.figure1,'UserData');
    v.contours{1} = v.contours{v.contourSelected};
    v.masks{1} = v.masks{v.contourSelected};
    %set(handles.figure1,'UserData',v);
    MasksChanged(handles);
    showImage(handles);

% --------------------------------------------------------------------
function Untitled_7_Callback(hObject, eventdata, handles)
global v
% hObject    handle to Untitled_7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    %v = get(handles.figure1,'UserData');
    v.contours{2} = v.contours{v.contourSelected};
    v.masks{2} = v.masks{v.contourSelected};
    %set(handles.figure1,'UserData',v);
    MasksChanged(handles);
    showImage(handles);

% --------------------------------------------------------------------
function Untitled_8_Callback(hObject, eventdata, handles)
global v
% hObject    handle to Untitled_8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    %v = get(handles.figure1,'UserData');
    v.contours{3} = v.contours{v.contourSelected};
    v.masks{3} = v.masks{v.contourSelected};
    %set(handles.figure1,'UserData',v);
    MasksChanged(handles);
    showImage(handles);


% --- Executes on mouse press over figure background.
function figure1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu8.
function popupmenu8_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu8 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu8
contents = cellstr(get(hObject,'String'));
set(handles.edit2,'String',contents{get(hObject,'Value')});

% --- Executes during object creation, after setting all properties.
function popupmenu8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
global v
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    saveRC(handles,0);
    %v = get(handles.figure1,'UserData');
    str = get(handles.edit2,'String');
    v.OldTimeStamp = now;
    if(abs(v.HitAKeyTimeStamp - v.OldTimeStamp) < 2e-5)
        return;
    end
    try
        set(handles.pushbutton10,'String','....');
        set(handles.edit2,'ForegroundColor',[0 0 0]);
        drawnow
        %set(handles.figure1,'UserData',v);
        copyfile(v.parfileName,'mpi2d.par');
        eval(str);
        loadCR(handles,1);
        showImage(handles);
        showCurves(handles);
        showShifts(handles);
    catch ME
        set(handles.edit2,'ForegroundColor',[153/256 0 0]);
        disp(ME.message);
    end
    set(handles.pushbutton10,'String','Eval');
    
    

% --- Executes on key press with focus on edit2 and none of its controls.
function edit2_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if(strcmp(eventdata.Key,'return'))
    pushbutton10_Callback(hObject, eventdata, handles);
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('Welcome to the Contour Regestratin Editor');
disp('Mouse')
disp('   Left drag = move frame');
disp('   Right drag = move Selected Contour');
disp('   Scroll = dilate/erode the Selected Contour');
disp('   Left drag(Bump Mode) = Bump in(blue circle while mouse is outside countour)');
disp('                               out(red circle while mouse is inside countour)');
disp('   Scroll(Bump Mode) = Enlarge/shrink the radius of the bump tool');
disp('Bump = Enable bump tool.  This is used to trim and extend contours');
disp('Pageup = go to previous frame');
disp('Pagedown = go to next frame');
disp('Space or p = play movie');
disp('Arrows = shift frame around');
disp('1 = Select the AIF');
disp('2 = Select the epicardium');
disp('3 = Select the endocardium');
disp('');
disp('mpi2d dropdown = preset list of commands to run mpi2d stages');
disp('mpi2d textbox = place to edit text to be evaluated');
disp('start angle = the insertion angle, or the most clockwise point in the RV');
disp('What are the k numbers floating around?   They are the k-trans values of the 2-compartment model');
disp('Below all the buttons in the upper right hand part of the screen are the shift vectors');
disp('Below the shift vectors are the tissue curves through time');


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
global showImageMutex bumpingMutex dragMutex keyPressMutex
keyPressMutex = 0;
dragMutex = 0;
showImageMutex = 0;
bumpingMutex = 0;
saveRC(handles,1);
delete(hObject);


% --- Executes on key press with focus on pushbutton10 and none of its controls.
function pushbutton10_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton10.
function pushbutton10_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
