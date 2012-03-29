function badfiles = FindSomething(directory)
global flaggedFiles
if(~exist('flaggedFiles'))
    flaggedFiles = {};
    flaggedFiles.badfiles = {};
end
files = dir(directory);
directoryCount = 1;
directories = {};
for i=3:length(files)
    if files(i).isdir
        directories(directoryCount) = {files(i).name};
        %disp([ char(directories(directoryCount))]);
        FindSomething([directory '/' char(directories(directoryCount))]);
        %disp(['Out  ' char(directories(directoryCount))]);
        directoryCount = directoryCount+1;
        
    end
end


files = dir([directory '/*.par']);
for i=1:length(files)
    fid = fopen([directory '/' char(files(i).name)]);
    
    tline = fgetl(fid);
    while ischar(tline)
        A = sscanf(tline,'%f %f');
        magnitude = sqrt(A(1)*A(1) + A(2)*A(2));
        if(magnitude > 10)
            %disp([num2str(magnitude) ' : ' directory '/' char(files(i).name)]);
            flaggedFiles.badfiles(end+1) = {[ directory '/' char(files(i).name)]};
            break;
        end
        tline = fgetl(fid);
    end
    fclose(fid);
end
badfiles = flaggedFiles.badfiles;