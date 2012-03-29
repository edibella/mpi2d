function filesFound = gatherFilesRecusivly()
    files = dir('*');
    filesFound = {};
    for i=3:length(files)
        if(~files(i).isdir)
            filesFound = vertcat(filesFound,files(i).name);
        else
            mypwd = pwd;
            cd(files(i).name);
            dug = gatherFilesRecusivly();
            for otherFilesi=1:length(dug)
                filesFound = vertcat(filesFound,[files(i).name '/' dug{otherFilesi}]);
            end
            cd(mypwd);
        end
    end
end

