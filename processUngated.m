function processUngated(raid1path, seriesName, first_series, numSlices, preBolusFlag, saveflag);
% to Process ungated perfusion.  9/19/11
% SET below parameters
% AFTER it is done, run 3.11 for all and fix, then mpi2d stg=[0.23 3.2 3.4 4.1
% 4.2]

% NOTE: Two functions put below in this file. Used to be separate scripts.


% reads in dicom like in 1.2, or
% read in cinemri, just recon with no header?
% find rv and lv centers.. 
% run self-gating

if preBolusFlag==1
    isub_range=[90 91]
else
    isub_range=[92 93]
end




for islice=1:numSlices  % chnage also dir name below!!
iseries=first_series+(islice-1)*1000;

seriesName_tmp=[seriesName,num2str(iseries),')']

numFrames=call_autoROI(raid1path, seriesName_tmp, islice, iseries, saveflag);

% cd([raid1path,'/Processing/']);
% 
% call_sys_dias_modifyParFile(raid1path, islice, iseries, isub_range, numFrames);
% 
% % next part is from Auto_all_ForDualBolus
%     
%     for isub=isub_range  % low dose sys, low dose dias, high sys, high dias
%       %  if isub==90 || isub==92
%       parfile=['series',int2str(iseries+isub),'.slice',int2str(islice),'.par'] 
%       if(exist('mpi2d.par')) delete('mpi2d.par'); end
%       copyfile(parfile,'mpi2d.par');
%     
%       disp(['Processing ' parfile]);
%     %2.6 cross-correlation registration with gaussian centering
%         %if this registration screws up, run mpi2d stg=2.11 to apply zero
%         %shifts then do 3.11 to do manual registration
%     %3.12 automatic segmentation
%     %3.2 extract curves
%     %3.4 convert curves to deltaSIcurves
%     %4.1 compute 2-compartment model params
%     %6.1 display k-trans
%       mpi2d stg=[2.8 3.12 3.2 3.4 4.1] myothreshold=.5
%     end


end





function call_sys_dias_modifyParFile(raid1path, islice, iseries, isub_range, numFrames)

templateParFile=['series',int2str(iseries),'.slice',int2str(islice),'.par'] 
    
%peakloc_complement=setdiff(1:max(peakloc),peakloc)
raid1path=[raid1path, '/ReconData/mat_files/'];
if ~exist(raid1path, 'dir')
    mkdir(raid1path)
end

flag_found_framesToSelect=0;
    
for isub=isub_range  % low dose sys, low dose dias, high sys, high dias
   if isub==90 || isub==92
       load([raid1path,'listofselectedFrames_roi_sum_dd_Nearsystolic_series',int2str(iseries),'slice',int2str(islice),'.mat'])
   else
       load([raid1path,'listofselectedFrames_roi_sum_dd_Neardiastolic_series',int2str(iseries),'slice',int2str(islice),'.mat'])
   end
   
   outfile=['series',int2str(iseries+isub),'.slice',int2str(islice),'.par']  % 91 is low dose diastole
   fidout=fopen(outfile,'w')
 
   fidin = fopen(templateParFile);
    if(fidin < 0)
        disp('I could not find the template file. Create with stg=1.2. It should be here:');
        disp(templateParFile);
    end
    while 1
        tline = fgetl(fidin);
        if ~ischar(tline)
            break;
        end
        tline
        %for each intervention check if this line has what we want
        %  if so write our own line
        if(~isempty(strfind(tline,'framesToSelect')))
            flag_found_framesToSelect=1;
            fprintf(fidout, 'framesToSelect=[');
            fprintf(fidout, '%d ', peakloc);
            fprintf(fidout, ']\n');
%         if(~isempty(strfind(tline,'ranget')))
%             fprintf(fidout, 'ranget=[');
%             fprintf(fidout, '%d ', peakloc);
%             fprintf(fidout, ']\n');
            
%         elseif(~isempty(strfind(tline,'ranget')))
%             fprintf(fidout, 'ranget=[1:110]\n');   % a guess at low dose range
        elseif(~isempty(strfind(tline,'ranget')))
            if numFrames < 250          % assume dual bolus if > 250
                fprintf(fidout, ['ranget=[7:',int2str(numFrames),']\n']);
            else
                if isub==90 || isub==91       
                   fprintf(fidout, 'ranget=[7:120]\n');  
                else
                   fprintf(fidout, ['ranget=[110:',int2str(numFrames),']\n']);  
                end   
            end
        elseif(~isempty(strfind(tline,'scaleAIF')))
            if isub==90 || isub==91
                fprintf(fidout, 'scaleAIF=1\n');  
            else
                %fprintf(fidout, 'scaleAIF=10\n'); 
                fprintf(fidout, 'scaleAIF=1\n'); % going to stick with own AIF first for now 10/6/11
            end           
        elseif(~isempty(strfind(tline,'studyNum')))
            fprintf(fidout, 'studyNum=%d\n',iseries+isub);   % a guess at low dose range
        elseif(~isempty(strfind(tline,'seriesNumAIF')))
            if isub==90 || isub==91
                fprintf(fidout, 'seriesNumAIF=%d\n', iseries+isub);
            else
                %fprintf(fidout, 'seriesNumAIF=%d\n', iseries+isub-2); 
                fprintf(fidout, 'seriesNumAIF=%d\n', iseries+isub);  % going to stick with own AIF first for now 10/6/11
            end
        else
            fprintf(fidout, [tline '\n'] );
        end
    end
    if flag_found_framesToSelect==0
        fprintf(fidout, 'framesToSelect=[');
        fprintf(fidout, '%d ', peakloc);
        fprintf(fidout, ']\n');
    end 
    fclose(fidin);
end

return



function numFrames=call_autoROI(raid1path, seriesName, islice, iseries, saveflag)
% call_autoROI, writes out nearSys and nearDias

dirname=strcat(raid1path,'ReconData/',seriesName,'/')

files=dir(strcat(dirname,'*.dcm'))
            for bb=1:length(files)
                tmpp=dicomread( strcat(dirname,files(bb).name));
                %tmpp=flipud(tmpp);
                %imgrr(:,:,bb)=rot90(tmpp);
                soss(:,:,bb)=double(tmpp);
            end
            

fullcinemri1=soss(:,:,100:round(bb/2));  % skip pre-bolus part, assuming it is less reliable?  Actually might be better! Should try. 
fullcinemri1=soss;

  [RV,LV] = FindLVRV(fullcinemri1,1);
  
  % just center on RV: 
  
  sx=size(fullcinemri1,1)
  sy=size(fullcinemri1,2)
  
            windowWidth = sx/3.3;
            windowHeight = sy/3.3;
            rangex = round(RV(1) - windowWidth/2):round(RV(1) + windowWidth/2);
            rangey = round(RV(2) - windowHeight/2):round(RV(2) + windowHeight/2);
            rangey = round(RV(2) - windowHeight/2):round(RV(2) + windowHeight);
            
            sz_adjust=1.5
            rangex2 = round(RV(1) - sz_adjust*windowWidth/2):round(RV(1) + sz_adjust*windowWidth/2);
            if round(RV(1) - sz_adjust*windowHeight/2) < 5
               rangex2=5:round(RV(1) + sz_adjust*windowHeight);
            end
            rangey2 = round(RV(2) - sz_adjust*windowHeight/2):round(RV(2) + sz_adjust*windowHeight);
            if round(RV(2) - sz_adjust*windowHeight/2) < 5
               rangey2=5:round(RV(2) + sz_adjust*windowHeight);
            end
            RV(1) = max(RV(1),min(rangex)+1);
            RV(2) = max(RV(2),min(rangey)+1);
            
            
  text(RV(2)+10, RV(1)-15,'RV','Color',[1 1 1]);
            plot([min(rangey) max(rangey) max(rangey) min(rangey) min(rangey)],[min(rangex) min(rangex) max(rangex) max(rangex) min(rangex)]);
            RV = RV - [min(rangex) min(rangey)]; %#ok<NASGU>
            LV = LV - [min(rangex) min(rangey)]; %#ok<NASGU>
            
         %keyboard   
            
numFrames=size(soss,3);
for ii=1:numFrames 
    ii
    
    currentImg=soss(rangex2,rangey2,ii);   %changed from rangex, EVRD 1/18/12, for P080411
   
    dd(ii)=sum(sum(currentImg));
   
    
end
    

 
dd_diff=diff(-dd);
peakloc=[];
for ii=1:length(dd_diff)-1
    if ( sign(dd_diff(ii))>sign(dd_diff(ii+1)) )
        peakloc=[peakloc ii+1];
    end
    
end
[ppvals, peaklocs2]=findpeaks(-dd, 'minpeakdistance',2)
selected_sos=soss(:,:,peaklocs2);
figNum=33;
figure(figNum); clf;
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',14)
dd=dd./max(dd);
lw=1.6
plot(dd,'k','linewidth',lw); hold on
plot(peaklocs2,dd(peaklocs2),'bo','linewidth',lw); hold on
%plot(locs,cc(locs),'mo','linewidth',lw); hold on
xlabel('Frame number')
legend('Sum in region ', 'Correlation in region ')
%whitebg('k')
%whitebg(figNum,[0.3 .3 0.3])


%title('dd, from sum ')
figure; 

imagesc(selected_sos(rangex,rangey,round(size(selected_sos,3)/2)))
figure;
%popd(selected_sos(rangex2,rangey2,:),' -dd minpeak2 from sum roi ')

% -dd seems more systolic
peakloc=peaklocs2   % fewer frames, otherwise the same. 
fileout=strcat(raid1path,'ReconData/mat_files/sum_roi_sum_minusdd_sys_series',int2str(iseries),'slice',int2str(islice),'.mat')
fileout=strcat(raid1path,'ReconData/mat_files/listofselectedFrames_sum_roi_minusdd_sys_series',int2str(iseries),'slice',int2str(islice),'.mat')
% save(fileout,'selected_sos')
if saveflag
    save(fileout, 'peakloc')
end



%adding 5/26/11, to write out the half closest to systole, then half
%closest to diastole.. 
systole_lowpeaks=peakloc;
[ppvals, peaklocsTmp]=findpeaks(dd, 'minpeakdistance',2)
diastole_highpeaks=peaklocsTmp;
new_diastole=[]; new_systole=[];
for ii=1:length(dd)
    % closest peak is:
           
    [minval_systole, indexsys]=min(abs(ii-systole_lowpeaks))
    [minval_diastole, indexdias]=min(abs(ii-diastole_highpeaks))
    
    if (minval_systole==0 || minval_diastole==0)
        if minval_systole==0
            new_systole=[new_systole ii];
        else
            new_diastole=[new_diastole ii];
        end
    else
        if ( abs(dd(ii)-dd(systole_lowpeaks(indexsys))) >= abs(dd(ii)-dd(diastole_highpeaks(indexdias))) )
            new_diastole=[new_diastole ii];
        else
            new_systole=[new_systole ii];
        end
    end
    
end

figure;
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',14)
plot(dd)
hold on
plot(new_systole,dd(new_systole),'ko','linewidth',lw)
plot(new_diastole,dd(new_diastole),'ro','linewidth',lw)
xlabel('Frame number')
ylabel('Arbitrary units')
legend('Sum in region ', 'Closest to systole', 'Closest to diastole')
selected_sos=soss(:,:,new_systole);
% figure; 
% popd(selected_sos(rangex2,rangey2,:),' new systole from sum roi ')
selected_sos=soss(:,:,new_diastole);
% figure; 
% popd(selected_sos(rangex2,rangey2,:),' new diastole from sum roi ')

fileout=strcat(raid1path,'ReconData/mat_files/listofselectedFrames_roi_sum_dd_Nearsystolic_series',int2str(iseries),'slice',int2str(islice),'.mat')
if saveflag
    peakloc=new_systole;
    save(fileout, 'peakloc')
end
fileout=strcat(raid1path,'ReconData/mat_files/listofselectedFrames_roi_sum_dd_Neardiastolic_series',int2str(iseries),'slice',int2str(islice),'.mat')
if saveflag
    peakloc=new_diastole;
    save(fileout, 'peakloc')
end

return


