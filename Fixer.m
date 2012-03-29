function varargout = Fixer(varargin)
% FIXER M-file for Fixer.fig
%      FIXER, by itself, creates a new FIXER or raises the existing
%      singleton*.
%
%      H = FIXER returns the handle to a new FIXER or the handle to
%      the existing singleton*.
%
%      FIXER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIXER.M with the given input arguments.
%
%      FIXER('Property','Value',...) creates a new FIXER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Fixer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Fixer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Fixer

% Last Modified by GUIDE v2.5 25-May-2010 08:17:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Fixer_OpeningFcn, ...
                   'gui_OutputFcn',  @Fixer_OutputFcn, ...
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


% --- Executes just before Fixer is made visible.
function Fixer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Fixer (see VARARGIN)

% Choose default command line output for Fixer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Fixer wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global v
clc
% if(exist('/v/raid1/bmatthew/mydata.mat'))
%     load('/v/raid1/bmatthew/mydata.mat');
%     v.patients = temp;
%     madeChanges = 0;
% else
%     clear temp
%     temp.shifts = [0];
%     v.patients = containers.Map({''},{temp});
%     remove(v.patients,'');
% end
root = '/v/raid1/bmatthew/MRIdata/Cardiac/Verio';
cd(root);
patients = dir([root '/*']);
for i=8%:length(patients)
    
    cd(root);
    parfiles = dir([patients(i).name '/Processing/*.par']);
    clear global v
    clear temp
    temp.shifts = [0];
    v.patients = containers.Map({''},{temp});
    remove(v.patients,'');
    
    for j=1:length(parfiles)
        if(isKey(v.patients,[patients(i).name ' ' parfiles(j).name]))
            continue;
        end
        if(strcmp(parfiles(j).name,'mpi2d.par') || ~isempty(strfind(parfiles(j).name,'rad')))
            continue;
        end
        disp('---------------');
        clear rangex rangey ranget infile studyNum sliceNum temp cinemri1
        ParFileName = parfiles(j).name;
        cd([patients(i).name '/Processing/']);
        ReadPar
        cd(root);
        temp.ParFileName = ParFileName;
        temp.rangex = rangex;
        temp.rangey = rangey;
        temp.ranget = ranget;
        temp.series = studyNum;
        temp.slice = sliceNum;
        temp.infile = infile;
        temp.shiftfileName = [patients(i).name '/Processing/Output/shiftsMAN.study' num2str(studyNum) '.slice' num2str(sliceNum) '.txt'];
        
        if(~exist( temp.shiftfileName))
            disp(['No shift file found for ' patients(i).name ' ' parfiles(j).name '. Putting zero shifts']);
            temp.shifts = [0 0];
        else
            temp.shifts = load( temp.shiftfileName);
        end
        
        if(~exist(temp.infile))
            disp(['Could not find infile for ' patients(i).name ' ' parfiles(j).name ': ' temp.infile]);
            continue;
        else
            if(~isempty(strfind(temp.infile,'.mat')))
                disp(['This refers to a mat file. Skipping. ' patients(i).name ' ' parfiles(j).name ': ' temp.infile]);
                continue;
            end
            temp.dicomfiles = dir([temp.infile '*.dcm']);
            if(length(temp.dicomfiles) < 40)
                disp([patients(i).name ' ' parfiles(j).name '  Less than 40 frames.  Skipping']);
                continue;
            else 
                disp([patients(i).name ' ' parfiles(j).name]);
            end
            for k=1:length(temp.dicomfiles)
                Header = dicominfo([temp.infile temp.dicomfiles(k).name]);
                temp.instance(k) = Header.InstanceNumber;
                temp.time(k) = str2num(Header.AcquisitionTime);
                temp.cinemri1original(:,:,k) = dicomread([temp.infile temp.dicomfiles(k).name]);
            end
            [temp.instance IX] = sort(temp.instance);
            temp.cinemri1original = temp.cinemri1original(:,:,IX);
            disp('Upsampling');
            temp.originalUpsampled = mpi_upsampleImages(temp.cinemri1original);
            temp.proposedLength = length(temp.shifts);
            tempshifts = zeros(size(temp.cinemri1original,3),2);
            for k=1:size(temp.shifts,1)
                tempshifts(k,:) = temp.shifts(k,:);
            end
            temp.shifts = tempshifts;
            for t=1:size(temp.originalUpsampled,3)
                frame = circshift(temp.originalUpsampled(:,:,t),temp.shifts(t,:));
                if(2*max(temp.rangex) > size(temp.originalUpsampled,1))
                    temp.originalExtracted(:,:,t) = frame(temp.rangex,temp.rangey);
                else
                    temp.originalExtracted(:,:,t) = frame((2*min(temp.rangex)):(2*min(size(frame,1),max(temp.rangex))),(2*min(temp.rangey)):(2*min(size(frame,2),max(temp.rangey))));
                end
            end
            temp.time = sort(temp.time);
            if(length(temp.shifts) > size(temp.cinemri1original,3))
                disp(['Found shift file that has more frames than original']);
                keyboard;
            end
            
        end
        
        temp.cinefilename = [patients(i).name '/Processing/Output/cinemri1.study' num2str(studyNum) '.slice' num2str(sliceNum) '.mat'];
        if(~exist(temp.cinefilename))
            disp(['No cine file found for ' patients(i).name ' ' parfiles(j).name '.  Using original']);
            temp.cinemri1 = temp.originalUpsampled;
        else
            load(temp.cinefilename);
            temp.cinemri1 = cinemri1;
            temp.proposedLength = max(temp.proposedLength,size(cinemri1,3));
        end
        
        temp.timefilename = [patients(i).name '/Processing/timeStampSer' num2str(studyNum) '.mat'];
        
        temp.bloodfilename = [patients(i).name '/Processing/Output/blood_polyCoords.study' num2str(studyNum) '.slice' num2str(sliceNum)];
        if(~exist(temp.bloodfilename))
            disp(['No blood file found for ' patients(i).name ' ' parfiles(j).name '.']);
            temp.useBlood = 0;
        else
            temp.useBlood = 1;
            temp.bloodcoords = load(temp.bloodfilename);
        end
        
        temp.epifilename = [patients(i).name '/Processing/Output/epi_polyCoords.study' num2str(studyNum) '.slice' num2str(sliceNum)];
        if(~exist(temp.epifilename))
            disp(['No epi file found for ' patients(i).name ' ' parfiles(j).name '.']);
            temp.useepi = 0;
        else
            temp.useepi = 1;
            temp.epicoords = load(temp.epifilename);
        end
        
        temp.endofilename = [patients(i).name '/Processing/Output/endo_polyCoords.study' num2str(studyNum) '.slice' num2str(sliceNum)];
        if(~exist(temp.endofilename))
            disp(['No endo file found for ' patients(i).name ' ' parfiles(j).name '.']);
            temp.useendo = 0;
        else
            temp.useendo = 1;
            temp.endocoords = load(temp.endofilename);
        end
        temp.t = 1;
        madeChanges = 1;
        v.patients([patients(i).name ' ' parfiles(j).name]) = temp;
        
    end
    if(madeChanges)
        disp('Saving');
        temp = v.patients;
        save(['/v/raid1/bmatthew/mydata' num2str(i) '.mat'],'temp');
    end
end
load('mydata8.mat');
global v
v.patients = temp;
v.handles = handles;
mylist = keys(v.patients);
v.current = mylist{1};
setCurrent(v.current);
set(handles.listbox1,'String',mylist);

%----------------------------
function setCurrent(whichOne)
global v
v.current = whichOne;
data = v.patients(v.current);

if(isfield(v,'bloodoriginal'))
    v = rmfield(v,'bloodoriginal');
end
if(isfield(v,'blood'))
    v = rmfield(v,'blood');
end

if(isfield(v,'curvesoriginal'))
    v = rmfield(v,'curvesoriginal');
end
if(isfield(v,'curves'))
    v = rmfield(v,'curves');
end

%----------------------------
axes(v.handles.axes1);
imagesc(data.cinemri1(:,:,data.t)), colormap gray
hold on
if(data.useBlood)
    plot(data.bloodcoords(:,1),data.bloodcoords(:,2));
end
if(data.useepi)
    plot(data.epicoords(:,1),data.epicoords(:,2));
end
if(data.useendo)
    plot(data.endocoords(:,1),data.endocoords(:,2));
end
hold off


%----------------------------
axes(v.handles.axes2);
imagesc(data.originalExtracted(:,:,data.t)), colormap gray
hold on
if(data.useBlood)
    plot(data.bloodcoords(:,1),data.bloodcoords(:,2));
end
if(data.useepi)
    plot(data.epicoords(:,1),data.epicoords(:,2));
end
if(data.useendo)
    plot(data.endocoords(:,1),data.endocoords(:,2));
end
hold off

%----------------------------
axes(v.handles.axes3);
if(size(data.originalExtracted,1) == size(data.cinemri1(:,:,data.t),1) && size(data.originalExtracted,2) == size(data.cinemri1(:,:,data.t),2))
    imagesc(abs(data.originalExtracted(:,:,data.t)-data.cinemri1(:,:,data.t))), colormap gray
else
    [s1 s2 s3] = size(data.originalExtracted);
    imagesc(rand(s1,s2));
end
hold on
if(data.useBlood)
    plot(data.bloodcoords(:,1),data.bloodcoords(:,2));
end
if(data.useepi)
    plot(data.epicoords(:,1),data.epicoords(:,2));
end
if(data.useendo)
    plot(data.endocoords(:,1),data.endocoords(:,2));
end
hold off
%----------------------------
axes(v.handles.axes4)
v.LVMask = roipoly(data.cinemri1(:,:,data.t),data.bloodcoords(:,1),data.bloodcoords(:,2));
v.epiMask = roipoly(data.cinemri1(:,:,data.t),data.epicoords(:,1),data.epicoords(:,2));
v.endoMask = roipoly(data.cinemri1(:,:,data.t),data.endocoords(:,1),data.endocoords(:,2));
for t=1:size(data.cinemri1,3)
    temp = data.cinemri1(:,:,t);
    v.blood(t) = mean(temp(v.LVMask>0));
end
[v.myoX v.myoY] = find(v.endoMask - v.epiMask);
v.myo = zeros(size(v.endoMask));
[x y] = find(v.LVMask);
v.center = mean([x y]);
for i=1:length(v.myoX)
    angle = atan2(v.myoX(i)-v.center(1),v.myoY(i)-v.center(2));
    r = norm([v.myoX(i)-v.center(1) v.myoY(i)-v.center(2)]);
    %profileBins(round(),round()) = v.data(
    region = floor(((angle+pi)/(2*pi)) * 6)+1;
    v.myo(v.myoX(i),v.myoY(i)) = region;
end

v.curves = zeros(6,size(data.cinemri1,3));
%myocardium
for i=1:6
    mask = (v.myo == i);
    for t=1:size(data.cinemri1,3)
        temp = data.cinemri1(:,:,t);
        v.curves(i,t) = mean(temp( mask>0));
    end
    v.curves(i,:) = v.curves(i,:)-mean(v.curves(i,1:7));
end

v.blood(:) = v.blood(:)-mean(v.blood(1:7));
plot(v.blood);
hold on
for i=1:6
    plot(v.curves(i,:));
end
hold off
%---------------------------
axes(v.handles.axes5)
v.LVMask = roipoly(data.originalExtracted(:,:,data.t),data.bloodcoords(:,1),data.bloodcoords(:,2));
v.epiMask = roipoly(data.originalExtracted(:,:,data.t),data.epicoords(:,1),data.epicoords(:,2));
v.endoMask = roipoly(data.originalExtracted(:,:,data.t),data.endocoords(:,1),data.endocoords(:,2));
for t=1:size(data.cinemri1original,3)
    temp = data.originalExtracted(:,:,t);
    v.bloodoriginal(t) = mean(temp(v.LVMask>0));
end
[v.myoX v.myoY] = find(v.endoMask - v.epiMask);
v.myo = zeros(size(v.endoMask));
[x y] = find(v.LVMask);
v.center = mean([x y]);
for i=1:length(v.myoX)
    angle = atan2(v.myoX(i)-v.center(1),v.myoY(i)-v.center(2));
    r = norm([v.myoX(i)-v.center(1) v.myoY(i)-v.center(2)]);
    %profileBins(round(),round()) = v.data(
    region = floor(((angle+pi)/(2*pi)) * 6)+1;
    v.myo(v.myoX(i),v.myoY(i)) = region;
end

v.curvesoriginal = zeros(6,size(data.cinemri1original,3));
%myocardium
for i=1:6
    mask = (v.myo == i);
    for t=1:size(data.cinemri1original,3)
        temp = data.originalExtracted(:,:,t);
        v.curvesoriginal(i,t) = mean(temp( mask>0));
    end
    v.curvesoriginal(i,:) = v.curvesoriginal(i,:)-mean(v.curvesoriginal(i,1:7));
end

v.bloodoriginal(:) = v.bloodoriginal(:)-mean(v.bloodoriginal(1:7));
plot(v.blood);
hold on
for i=1:6
    plot(v.curvesoriginal(i,:));
end
hold off

%-----------------------------
axes(v.handles.axes6)
tmax=min(length(v.curvesoriginal),length(v.curves));
for i=1:6
    plot(abs(v.curvesoriginal(i,1:tmax) - v.curves(i,1:tmax)));
    hold on;
end
hold off;

%----------------------------
axes(v.handles.axes7);
clear temp
for t=1:length(data.shifts)
    temp(t) = norm(data.shifts(t,:));
end
plot(temp);

% --- Outputs from this function are returned to the command line.
function varargout = Fixer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global v
data = v.patients(v.current);
data.shifts = circshift(data.shifts,1);
for t=1:size(data.originalUpsampled,3)
    frame = circshift(data.originalUpsampled(:,:,t),data.shifts(t,:));
    data.originalExtracted(:,:,t) = frame((2*min(data.rangex)):(2*min(size(frame,1),max(data.rangex))),(2*min(data.rangey)):(2*min(size(frame,2),max(data.rangey))));
end
v.patients(v.current) = data;
 setCurrent(v.current);

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global v
data = v.patients(v.current);
data.shifts = circshift(data.shifts,-1);
for t=1:size(data.originalUpsampled,3)
    frame = circshift(data.originalUpsampled(:,:,t),data.shifts(t,:));
    data.originalExtracted(:,:,t) = frame((2*min(data.rangex)):(2*min(size(frame,1),max(data.rangex))),(2*min(data.rangey)):(2*min(size(frame,2),max(data.rangey))));
end
v.patients(v.current) = data;
 setCurrent(v.current);

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
contents = get(hObject,'String');
setCurrent(contents{get(hObject,'Value')});


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global v
data = v.patients(v.current);
temp = data.shifts(:,1);
data.shifts(:,1) = data.shifts(:,2);
data.shifts(:,2) = temp;
for t=1:size(data.originalUpsampled,3)
    frame = circshift(data.originalUpsampled(:,:,t),data.shifts(t,:));
    data.originalExtracted(:,:,t) = frame((2*min(data.rangex)):(2*min(size(frame,1),max(data.rangex))),(2*min(data.rangey)):(2*min(size(frame,2),max(data.rangey))));
end
v.patients(v.current) = data;
setCurrent(v.current);

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over listbox1.
function listbox1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global v
ratio = (get(hObject,'Value')-get(hObject,'Min'))/get(hObject,'Max');
data = v.patients(v.current);
data.t = ceil(max(size(data.cinemri1,3),size(data.originalUpsampled,3))*ratio);

%----------------------------
axes(v.handles.axes1);
imagesc(data.cinemri1(:,:,data.t)), colormap gray
hold on
if(data.useBlood)
    plot(data.bloodcoords(:,1),data.bloodcoords(:,2));
end
if(data.useepi)
    plot(data.epicoords(:,1),data.epicoords(:,2));
end
if(data.useendo)
    plot(data.endocoords(:,1),data.endocoords(:,2));
end
hold off


%----------------------------
axes(v.handles.axes2);
imagesc(data.originalExtracted(:,:,data.t)), colormap gray
hold on
if(data.useBlood)
    plot(data.bloodcoords(:,1),data.bloodcoords(:,2));
end
if(data.useepi)
    plot(data.epicoords(:,1),data.epicoords(:,2));
end
if(data.useendo)
    plot(data.endocoords(:,1),data.endocoords(:,2));
end
hold off

%----------------------------
axes(v.handles.axes3);
if(size(data.originalExtracted,1) == size(data.cinemri1(:,:,data.t),1) && size(data.originalExtracted,2) == size(data.cinemri1(:,:,data.t),2))
    imagesc(abs(data.originalExtracted(:,:,data.t)-data.cinemri1(:,:,data.t))), colormap gray
else
    [s1 s2 s3] = size(data.originalExtracted);
    imagesc(rand(s1,s2));
end
hold on
if(data.useBlood)
    plot(data.bloodcoords(:,1),data.bloodcoords(:,2));
end
if(data.useepi)
    plot(data.epicoords(:,1),data.epicoords(:,2));
end
if(data.useendo)
    plot(data.endocoords(:,1),data.endocoords(:,2));
end
hold off

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on key press with focus on slider1 and none of its controls.
function slider1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global v
data = v.patients(v.current);
data.originalUpsampled= circshift(data.originalUpsampled,[0 0 1]);
for t=1:size(data.originalUpsampled,3)
    frame = circshift(data.originalUpsampled(:,:,t),data.shifts(t,:));
    data.originalExtracted(:,:,t) = frame((2*min(data.rangex)):(2*min(size(frame,1),max(data.rangex))),(2*min(data.rangey)):(2*min(size(frame,2),max(data.rangey))));
end
v.patients(v.current) = data;
setCurrent(v.current);

% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global v
data = v.patients(v.current);
data.originalUpsampled= circshift(data.originalUpsampled,[0 0 -1]);
for t=1:size(data.originalUpsampled,3)
    frame = circshift(data.originalUpsampled(:,:,t),data.shifts(t,:));
    data.originalExtracted(:,:,t) = frame((2*min(data.rangex)):(2*min(size(frame,1),max(data.rangex))),(2*min(data.rangey)):(2*min(size(frame,2),max(data.rangey))));
end
v.patients(v.current) = data;
setCurrent(v.current);
