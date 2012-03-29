clear 
clc
top = '/v/raid1/bmatthew/MRIdata/Cardiac/Trio';
patientFolders = dir([top '/*']);
for i=3:length(patientFolders)
    if(~patientFolders(i).isdir) continue; end;
    cd([top '/' patientFolders(i).name '/Processing']);
    parfiles = dir('*.par');
    disp(patientFolders(i).name);
    for j=1:length(parfiles)
        if(~isempty(strfind(parfiles(j).name,'rad'))) continue;  end;
        disp(parfiles(j).name);
        ParFileName = parfiles(j).name;
        ReadPar
        if(~isempty(strfind(infile,'.mat')))
            disp('This par file refers to a mat file');
            continue;
        end
        dicomFiles = dir([infile '*.dcm']);
        if(length(dicomFiles) == 0)
            disp('Found no dicom files');
        end
        if(length(dicomFiles)-1 ~= length(ranget))
            disp('ranget does not match dicom files');
        end
        shiftfilename = ['Output/shiftsMAN.study' num2str(studyNum) '.slice' num2str(sliceNum) '.txt'];
        if(~exist(shiftfilename))
            disp('Could not find shift file');
            continue;
        end
        shifts = load(shiftfilename);
        if(length(shifts) ~= length(ranget))
            if(abs(length(shifts) - length(ranget)) ==1)
               disp('Fixing shifts length by 1');
               %shifts = shifts(2:end,:);
               copyfile([shiftfilename '.backup'],shiftfilename);
               %fid = fopen(shiftfilename);
               %for t=1:length(shifts)
               %    fprintf(fid,'%f %f\n',shifts(t,1),shifts(t,2));
               %end
               %fclose(fid);
            else
                disp(['Difference = ' num2str(abs(length(shifts) - length(ranget)))]);
            end
            
            disp('shifts size does not match ranget');
        end
        load(['Output/cinemri1.study' num2str(studyNum) '.slice' num2str(sliceNum) '.mat']);
        if(size(cinemri1,3) ~= length(shifts))
            disp('Cinemri1 time length does not match shifts');
        end
        
    end
end