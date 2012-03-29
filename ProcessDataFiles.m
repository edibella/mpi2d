function ProcessDataFiles(infolder, outfolder, study,MID)
    if(infolder(end) ~= '/')
        infolder = [infolder '/'];
    end
    if(outfolder(end) ~= '/')
        outfolder = [outfolder '/'];
    end
    s = dir(infolder);
    n = length(s);
    for i=3:n
    	%values = sscanf(s.mat(i), '%*s%d%*s%d_slice%d_subset%d_%diter_%dcoils_F%d_T%d_S%d.mat', [1, inf, inf, inf, inf, inf, inf, inf, inf, inf])
        infile = s(i).name;
        if exist('MID')
            if(size(strfind(infile,['MID' int2str(MID)]))==0)
                continue;
            end
        end
        if(~(size(strfind(infile,'slice'))>0))
            continue;
        end
        s2 = split(s(i).name, '_');
        if size(s2,2) == 1      % was it not _ delimited?
            s2 = split(s(i).name,'.');
        end
        subset = -1; % default value if no subset is found
        slice = -1; % default value if no subset is found
        disp(['converting ' infile]);
        for j=1:size(s2,2)
            temp = char(s2(j));
            if(size(strfind(temp,'slice'))>0)
                slice = sscanf(temp,'slice%d',[1]);
            end
            if(size(strfind(temp,'x'))>0)
                dims = split(temp,'x');
                tmax = char(dims(3));
            end
            if exist('MID')
                if(size(strfind(temp,'subset'))>0)
                    subset = sscanf(temp,'subset%d',[1]);
                end
            end
        end
        if( slice < 0)
          disp('I could not extract the slice information from the following filename:');
          disp(infile);
        end
        makeParFile(infolder, infile, outfolder, '/v/raid1/bmatthew/Template.par', slice, subset, study,tmax)
    end
    %values = sscanf(mixed(k,:), '%*s %d %d %*s', [1, inf])