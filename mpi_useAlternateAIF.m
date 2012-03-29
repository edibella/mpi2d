function deltaSIcurves=mpi_useAlternateAIF(seriesNumAIF, sliceNumAIF, scaleAIF, sliceNum,seriesNum, outpath, numPreContrast_bld, numPreContrast_tiss, numSkip,framesToSelect, ranget, timeStampFile, timeStampFileAIF)

numCCs=scaleAIF;
delta_t=0.5
% if pre-bolus has fewer time frames just replace time frames up to that point...  7/5/05


%disp('WARNING!!!  THIS WILL OVERWRITE THE CURVES FILES TO PUT IN ALTERNATE AIF. SCALING MAY BE OFF!!!  Edit mpi2d if needed.')
%disp('Hit CNTRL-C to abort if you do not want curves overwritten! Easy to get back by re-running 3.2 ')
%studyNumAIF=input('Enter study or series number to get AIF from \n')
%sliceNumAIF=input('Enter slice number to get AIF from \n')

% read in alternate AIF 
   curvefilename=strcat(outpath,'deltaSIcurves.study',int2str(seriesNumAIF),'.slice',int2str(sliceNumAIF),'.mat')
   try load(curvefilename);
   catch disp('cant seem to find curves.series file to load, may not have run 3.2 yet')
   end
   aif_curve=deltaSIcurves(1,:)';
   [nRegsAIF, nTimesAIF]=size(deltaSIcurves);

 first_deltaSI=deltaSIcurves;
   clear deltaSIcurves

   curvefilename=strcat(outpath,'deltaSIcurves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat');
   load(curvefilename);
   second_deltaSI=deltaSIcurves;

   si_curves = deltaSIcurves;
   [nRegs, nTimes]=size(deltaSIcurves);
   nRegs=nRegs-1;
   bldcurve=si_curves(1,:)';

   % scale last 5 points of AIF to match?
   init_aif=sum(aif_curve(numSkip+1:numPreContrast_bld+numSkip))/numPreContrast_bld;
   init_bld=sum(bldcurve(numSkip+1:numPreContrast_bld+numSkip))/numPreContrast_bld;
%   bldcurve=aif_curve*10*init_bld/init_aif;





%%% To scale the AIF from this slice to another slice(since above the mean
%%% pre-contrast frames are near 0, I'll use the original curves files to get the pre-contrast scale factor....
%%% NP 04/27/07 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp1=strcat(outpath,'curves.study',int2str(seriesNumAIF),'.slice',int2str(sliceNumAIF),'.mat');
load(tmp1);
tmpAIFcurve=curves(1,:)';
tmp2=strcat(outpath,'curves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat');
load(tmp2);
tmpBLDcurve=curves(1,:)';
init_aif=sum(tmpAIFcurve(numSkip+1:numPreContrast_bld+numSkip))/numPreContrast_bld;
init_bld=sum(tmpBLDcurve(numSkip+1:numPreContrast_bld+numSkip))/numPreContrast_bld;
%%% to do this, scale_to_AIF should be set to 1 below!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% read in ranget and framesToSelect from   - should do numSkip too and
% framesToSkip
flag_framesToSelect=0

ParFileName=strcat('series',int2str(seriesNumAIF),'.slice',int2str(sliceNumAIF),'.par');
    fidin = fopen(ParFileName);
    if(fidin < 0)
        disp('I could not find the AIF par file. It should be here:');
        disp(ParFileName);
    end
    while 1
        tline = fgetl(fidin);
        if ~ischar(tline)
            break;
        end
        if(~isempty(strfind(tline,'ranget')))
            ranget_AIF=str2num(tline(strfind(tline,'=')+1:end))
        end
        if(~isempty(strfind(tline,'framesToSelect')))
            framesToSelect_AIF=str2num(tline(strfind(tline,'=')+1:end))
            flag_framesToSelect=1;% if there is no framesToSelect, create it from framesToSkip.
            % Easier - just use ranget... 
        end
    end
    if flag_framesToSelect==0;
        framesToSelect_AIF=ranget_AIF;
    end
            
% resample to uniform time points:
% first do low dose AIF, if present:
flagTimeStamps=0;
if exist('timeStampFileAIF')
   flagTimeStamps=1;
   load(timeStampFileAIF);   % assumes gives variable timeStamp
   if(size(timeStamp,1) == 1)
       %the timestamp was not saved in the format mpi2d expects
       timeStamp = timeStamp';
       save( timeStampFile,'timeStamp');
   end
    counter=1;
        for j=framesToSelect_AIF  %ranget    % EVRD 4/26/11, to crop to ranget. Assume subset of shifts. 
            if (j <=max(ranget_AIF)) && (j >=min(ranget_AIF))               
                selectedTimeStamps(counter,1) = timeStamp(j);
                counter=counter+1;
            end
        end
    
   aif_curve=interpTimeCurve(selectedTimeStamps, delta_t, aif_curve);
   %nTimes=length(bldcurve);

   
end


clear selectedTimeStamps
flagTimeStamps=0;
if exist('timeStampFile')
   flagTimeStamps=1;
   load(timeStampFile);   % assumes gives variable timeStamp
   if(size(timeStamp,1) == 1)
       %the timestamp was not saved in the format mpi2d expects
       timeStamp = timeStamp';
       save( timeStampFile,'timeStamp');
   end
    counter=1;
    if isempty(framesToSelect)
        framesToSelect=ranget;
    end
    
        for j=framesToSelect  %ranget    % EVRD 4/26/11, to crop to ranget. Assume subset of shifts. 
            if (j <=max(ranget)) && (j >=min(ranget))               
                selectedTimeStamps(counter,1) = timeStamp(j);
                counter=counter+1;
            end
        end
    
   bldcurve=interpTimeCurve(selectedTimeStamps, delta_t, bldcurve);
   sat_bldcurve=bldcurve;
   nTimes=length(bldcurve);

% now do different one for Tissues:  
   tisscurves = si_curves(2:nRegs+1,:);
   tisscurve=tisscurves';
   clear tisscurves;
   
   for ii=1:nRegs
       ii
       try
           tisscurves(ii,:)=interpTimeCurve(selectedTimeStamps, delta_t, tisscurve(:,ii));
        catch
            disp('Likely a problem with a tissue region having Nans. Take a look with 3.21 or 3.11, redraw ')
        end
   end
%   sat_bldcurve=interpTimeCurve(selectedTimeStamps, delta_t, sat_bldcurve);

% add sanity check to make sure AIF and tissue curves are same length, since they
% could come from different timeStamp files
%    if (length(bldcurve) > size(tisscurve,1))
%        bldcurve=bldcurve(1:size(tisscurve,1));  % tr
%        disp('Bld and tiss curves lenghts do not match due to timestamps...\n')
%        disp('Truncating bldcurve length to match tisscurve')
%        nTimes=size(tisscurve,1);
%    end
%    if (length(bldcurve) < size(tisscurve,1))
%        tisscurve=tisscurve(1:nTimes,:);
%        disp('Bld and tiss curves lenghts do not match due to timestamps...\n')
%        disp('Truncating tisscurve length to match bldcurve')
%    end
end  







   disp('USING ALTERNATE AIF')
scale_to_AIF=1; % set to 0 so scaling AIF directly (when volume matched) NP 041709 % this is set to 1 only if scaling AIF from another slice
add_inputsFlag=0;
summedAIF=0;


tmp_aif_curve=[aif_curve; zeros(length(sat_bldcurve)-length(aif_curve),1)];
[estGlobalDelay, tmp_aif_curve]=alignCurves(10*tmp_aif_curve, sat_bldcurve);  % shifts bldcurve to match sat_bldcurve, integer shifts. So can use single delay in fits.
% just use estGlobalDelay, not the shifted zero padded curve


if (scale_to_AIF==1) 
% don't scale up by 10, just scale resulting param. ktrans by 10...
  % bldcurve(1:nTimesAIF)=aif_curve*init_bld/init_aif;
  % disp('scaling aif by init_bld/init_aif NOT times 10')
  % init_bld/init_aif
   %bldcurve(1:length(aif_curve))=aif_curve*scaleAIF;
   try 
       bldcurve(-(estGlobalDelay+1):length(aif_curve)-(estGlobalDelay+2))=aif_curve*scaleAIF;
   catch
       bldcurve(1:length(aif_curve))=aif_curve*scaleAIF;
   end
   
else
   bldcurve(1:length(aif_curve))=aif_curve*numCCs;
   %%%%% if slice num doesn't match, mult by :  %%%%% bldcurve(1:nTimesAIF)=aif_curve*init_bld/init_aif; too
   if add_inputsFlag==1
      delta_t=0.6;
      disp('USING HARD CODED DELTA_T HERE!!!!')
      for k=1:numCCs-1
        t_delay=-k/6;   % cc/sec
        t_delay = t_delay/delta_t;   % convert from seconds to frames
        tmpbldcurve(1-fix(t_delay):nTimes)=bldcurve(1:nTimes+fix(t_delay));
        tmpbldcurve(1:1-fix(t_delay))=0;
        t_delay=t_delay-fix(t_delay);
        tmpbldcurve = interp1(1:nTimes,tmpbldcurve,1+t_delay:nTimes+t_delay,'splines',0);
        summedAIF=tmpbldcurve + summedAIF;
      end
      bldcurve(1:nTimesAIF)=summedAIF(1:nTimesAIF);
    end

end
%    tisscurves = si_curves(2:nRegs+1,:);
%    tisscurve=tisscurves';

   %%%% In both bldcurve lines below, only take range from 1:nTimes of low-
   %%%% dose AIF so lengths are the same.  This won't work if the length of
   %%%% the low-dose AIF is LESS than the new high-dose AIF   NP 10/27/06
deltaSIcurves =[ bldcurve(1:nTimes)'; tisscurves; sat_bldcurve'];  % so in same format as curves

showcurves([bldcurve(1:nTimes)'; tisscurves],'Alternate AIF shown');
return;



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


function res=delspace(str)

 in=find(str~=' ');
 if(~isempty(in)), res=str(in); end;

 
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

