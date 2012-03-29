function tissViewer(sliceToSee)
%% things to change
slicesToAllow = [1 2 3 4 5 6 7];% 3 4 5 6 7 8];

%% adjust ratios here
% run the program and adjust the ratios here by injection amount
injectioratio = [5 5 3.6/.7 1 2/3 2/3 2/3];
%injectioratio = [5 1 2/3 2/3 2/3];

injectioratio = [5 5 1 2/3 2/3 2/3]; %P041910
injectioratio = [1 1  1 1 1 1]; 

injectionsToShow = '[1 2 3 4]';
injectionsToShow = '1:length(injectionTypes)';

lineTypes = {'--*','--o','--s','--d','--^','--v','-->','--<'};


%% gather the curves based on the par files avaliable
mypwd = pwd;
if(isempty(strfind(mypwd,'Processing')))
    if(exist('Processing') == 7)
        cd('Processing');
    else
        disp('I cannot find the Processing folder.  Please naviagte to it');
        return;
    end
end
if(~isempty(strfind(mypwd,'Output')))
    cd('..');
end
parfiles = dir('*.par');
curveFileNames = {};
clear temp; temp.hello = 0;
mymap = containers.Map({''},{temp});
remove(mymap,'');

for p=1:length(parfiles)
    if(strcmp(parfiles(p).name,'mpi2d.par')) continue; end    %kill off mpi2d
    if(~isempty(strfind(parfiles(p).name,'rad'))) continue; end   %kill off subsets
    
    %extract the series and slice number
    A = sscanf(parfiles(p).name,'series%d.slice%d.par');
    if(isempty(A)) continue; end
    series = A(1); slice = A(2); %subset = A(3);
    if(isempty(find(slicesToAllow == slice)))
        continue;
    end
    
    %get at the injection type
    ParFileName = parfiles(p).name;
    ReadPar
    where = strfind(infile,'/');
    dicomFolder = infile((where(end-1)+1):(where(end)-1));
    parenStart = strfind(dicomFolder,'(');
    injectionType = dicomFolder(1:(parenStart-1));
    
    %find the delteSIcurves file
    template = strrep(strrep(parfiles(p).name,'series','study'),'par','mat');
    curvesCandidates = dir(['Output/deltaSIcurves*' template])
    
    %throw it in the map
    if(~isKey(mymap,injectionType))
        clear temp
        load(['Output/' curvesCandidates(1).name]);
        temp.curves = containers.Map({[0]},{deltaSIcurves});
        temp.series = containers.Map({[0]},{[0]});
        remove(temp.curves,[0]);
        mymap(injectionType) = temp;
    end
    data = mymap(injectionType);
    data.curvesFileName = curvesCandidates(1).name;
    load(['Output/' curvesCandidates(1).name]);
    data.curves(slice) = deltaSIcurves;
    data.series(slice) = series
end

disp('------------------');
%% display the curves
figure(44);
%set(gcf,'Color',[.23,0.44,0.34])
hold on
legendStrings = {};
%mycolormap = lines;    % change this to be whatever colormap you want the curves to divide between themselves
mycolormap = hsv;   %the curves are evenly spaced using the colormap given here
%mycolormap = gray;   %use this one if you want to print things out and kill the 'Color',[.23 ,.44...] thing above 
injectionTypes = keys(mymap);
for injectioni = eval(injectionsToShow)
    injection = injectionTypes{injectioni};
    disp(['Ratio(' injection ') = ' num2str(injectioratio(injectioni))])
    
    data = mymap(injection);
    slices = keys(data.curves);
%    seriesplusslices = keys(data.series)
    %for slicei = sliceToSee  %1:length(slices)
    for slicei = 1:length(slices)
        slice = slices{slicei};
        curves = data.curves(slice);
        tisscurve=mean(curves(2:end,:),1);
        AIF = curves(1,:);
        [maxc maxi] = max(AIF);
        offset = 22-maxi;
        AIF(1:max(1,offset)) = 0;
        AIF = circshift(AIF,[0 offset]);
        AIF(1:max(1,offset)) = 0;
        AIF = AIF*injectioratio(injectioni);
        tisscurve = tisscurve*injectioratio(injectioni);
        plot(tisscurve,lineTypes{slicei},'Color',1-mycolormap(1+floor((injectioni-1)*length(mycolormap)/length(injectionTypes)),:),'LineWidth',.5);
        legendStrings{length(legendStrings)+1} = [injection '(' num2str(data.series(slice)) ',' num2str(slice) ')'];
    end
end
hold off
xlabel('Time Frame');
ylabel('Intensity (delta SI)');
title({'AIF seperated by color(injection) and marker type(slice)'})
hold off
if(isempty(legendStrings)) 
    return; 
end
lineToExpress = ['legend(''' strrep(legendStrings{1},'_','\_') ''''];
for line=2:length(legendStrings)
    lineToExpress = [lineToExpress ',' '''' strrep(legendStrings{line},'_','\_') ''''];
end
lineToExpress = [lineToExpress ');'];
eval(lineToExpress);
end