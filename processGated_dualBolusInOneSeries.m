function processGated_dualBolusInOneSeries(raid1path, seriesName, first_series, numSlices, preBolusFlag);
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
    isub_range=[99]
else
    isub_range=[]
end




for islice=1:numSlices  % chnage also dir name below!!
iseries=first_series+(islice-1)*1000;

seriesName_tmp=[seriesName,num2str(iseries),')'];

dirname=strcat(raid1path,'ReconData/',seriesName_tmp,'/');

files=dir(strcat(dirname,'*.dcm'));
numFrames=length(files)
      
cd([raid1path,'/Processing/']);

call_makeParFiles_dualBolus(raid1path, islice, iseries, isub_range, numFrames);

% next part is from Auto_all_ForDualBolus
    
    for isub=isub_range  % low dose sys, low dose dias, high sys, high dias
      %  if isub==90 || isub==92
      parfile=['series',int2str(iseries+isub),'.slice',int2str(islice),'.par'] 
      if(exist('mpi2d.par')) delete('mpi2d.par'); end
      copyfile(parfile,'mpi2d.par');
    
      disp(['Processing ' parfile]);
    %2.6 cross-correlation registration with gaussian centering
        %if this registration screws up, run mpi2d stg=2.11 to apply zero
        %shifts then do 3.11 to do manual registration
    %3.12 automatic segmentation
    %3.2 extract curves
    %3.4 convert curves to deltaSIcurves
    %4.1 compute 2-compartment model params
    %6.1 display k-trans
      mpi2d stg=[2.8 3.12 3.2 3.4 4.1] myothreshold=.5
    end


end





function call_makeParFiles_dualBolus(raid1path, islice, iseries, isub_range, numFrames)

templateParFile=['series',int2str(iseries),'.slice',int2str(islice),'.par'] 
    
%peakloc_complement=setdiff(1:max(peakloc),peakloc)
raid1path=[raid1path, '/ReconData/mat_files/'];
if ~exist(raid1path, 'dir')
    mkdir(raid1path)
end

flag_found_framesToSelect=0;
    
for isub=isub_range  % low dose sys, low dose dias, high sys, high dias
   
   outfile=['series',int2str(iseries+isub),'.slice',int2str(islice),'.par']  % 99 is low dose
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

     
%         elseif(~isempty(strfind(tline,'ranget')))
%             fprintf(fidout, 'ranget=[1:110]\n');   % a guess at low dose range
        if(~isempty(strfind(tline,'ranget')))
            if numFrames < 150          % assume dual bolus if > 150
                fprintf(fidout, ['ranget=[7:',int2str(numFrames),']\n']);
            else
                if isub==99      
                   fprintf(fidout, 'ranget=[7:80]\n');  
                else
                   fprintf(fidout, ['ranget=[90:',int2str(numFrames),']\n']);  
                end   
            end
        elseif(~isempty(strfind(tline,'scaleAIF')))
            if isub==99 
                fprintf(fidout, 'scaleAIF=1\n');  
            else
                %fprintf(fidout, 'scaleAIF=10\n'); 
                fprintf(fidout, 'scaleAIF=1\n'); % going to stick with own AIF first for now 10/6/11
            end           
        elseif(~isempty(strfind(tline,'studyNum')))
            fprintf(fidout, 'studyNum=%d\n',iseries+isub);   % a guess at low dose range
        elseif(~isempty(strfind(tline,'seriesNumAIF')))
            if isub==99
                %fprintf(fidout, 'seriesNumAIF=%d\n', iseries - 1000*islice +isub);
                fprintf(fidout, 'seriesNumAIF=%d\n', iseries +isub);
            else
                fprintf(fidout, 'seriesNumAIF=%d\n', iseries+isub);  % going to stick with own AIF first for now 10/6/11
            end
        elseif(~isempty(strfind(tline,'sliceNumAIF')))
            if isub==99
                fprintf(fidout, 'sliceNumAIF=%d\n', islice);
            else
                %fprintf(fidout, 'seriesNumAIF=%d\n', iseries+isub-99); 
                fprintf(fidout, 'sliceNumAIF=%d\n', islice);  % going to stick with own AIF first for now 10/6/11
            end
        else
            fprintf(fidout, [tline '\n'] );
        end
    end
   
    fclose(fidin);
end

return


