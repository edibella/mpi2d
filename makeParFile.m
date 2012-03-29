function file = makeParFile(infolder, outfolder, templateParFile, slice, subset, series,tmax, rangex, rangey)
    %load([infolder infile]);
    if(~strcmp(outfolder(end),'/'))
        outfolder = [outfolder '/'];
    end
    if(~strcmp(infolder(end),'/'))
        infolder = [infolder '/'];
    end
    if(subset > 0)
        file = sprintf( 'rad.series%d.slice%d.subset%d.par',series,slice,subset);
    else
        file = sprintf('series%d.slice%d.par',series,slice);
    end
    [fid message] = fopen([outfolder file] , 'w');
    if(fid < 0)
        folders = split(outfolder,'/');
        current = ['/' char(folders(1))];
        makeDirectories = 0;
        for j=2:size(folders,2)
            if(isempty(char(folders(j)))) 
                break;
            end
            if(exist(current,'dir') ~= 7)
                makeDirectories = 1;
            end
            if(makeDirectories > 0)
                mkdir(current,char(folders(j)));
                if(exist(current,'dir') ~= 7)
                    disp(['I tried to make a folder named:' char(folders(j)) 'but could not']);
                    disp(['with the following message ' message]);
                    return;
                end
            end
            current = [current '/' char(folders(j))];
           
        end
        [fid message] = fopen([current '/' file] , 'w');
        if(fid < 0)
            disp('I tried to make the file but failed with this message:');
            disp(message);
            keyboard;
            return;
        end
        %disp('I could not open the output folder.  Please make it for me.  Pretty pretty please');
    end
    %tmax = size( imgrr,3);
    fidin = fopen(templateParFile);
    if(fidin < 0)
        disp('I could not find the template file. It should be here:');
        disp(templateParFile);
    end
    while 1
        tline = fgetl(fidin);
        if ~ischar(tline)
            break;
        end
        %for each intervention check if this line has what we want
        %  if so write our own line
        if(~isempty(strfind(tline,'seriesNumAIF')))
            fprintf(fid, 'seriesNumAIF=%d\n', series);
        elseif(~isempty(strfind(tline,'sliceNumAIF')))
            fprintf(fid, 'sliceNumAIF=%d\n', slice);
        elseif(~isempty(strfind(tline,'studyNum')))
            fprintf(fid, 'studyNum=%d\n', series);
        elseif(~isempty(strfind(tline,'sliceNum')))
            fprintf(fid, 'sliceNum=%d\n', slice);
        elseif(~isempty(strfind(tline,'infile=')))
            fprintf(fid, 'infile=%s\n', infolder);
        elseif(exist('rangex') && (~isempty(strfind(tline,'rangex='))))
            fprintf(fid, 'rangex=%d:%d\n', min(rangex),max(rangex));
        elseif(exist('rangey') && (~isempty(strfind(tline,'rangey='))))
            fprintf(fid, 'rangey=%d:%d\n', min(rangey),max(rangey));
        elseif(~isempty(strfind(tline,'timeStampFile')))
            fprintf(fid, 'timeStampFile=timeStampSer%d.mat\n',series );
        elseif(~isempty(strfind(tline,'ranget')))
            fprintf(fid, 'ranget=2:%d\n', tmax);
        elseif(~isempty(strfind(tline,'lastFrame')))
            fprintf(fid, 'lastFrame=%d\n', tmax-1);
        else 
            fprintf(fid, [tline '\n'] );
        end
    end
    fclose(fidin);fclose(fid);
end