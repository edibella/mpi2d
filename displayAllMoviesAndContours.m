function displayAllMoviesAndContours()
    parfiles = dir('*.par');
    outpath='Output/'  %EVRD 12/22/10
    previousInfileName = '';
    injections = containers.Map({''},{{''}}); remove(injections,'');
    for i=1:length(parfiles)
        if(~isempty(strfind(parfiles(i).name,'rad')) || strcmp(parfiles(i).name,'mpi2d.par'))
            continue;
        end
        %get at the series number and slice number
        ParFileName = parfiles(i).name %#ok<NASGU>
        ReadPar;
        
        %group them into injections
        genericInfilename = infile(1:strfind(infile,'('))
        if(~strcmp(genericInfilename,previousInfileName))
            previousInfileName = genericInfilename;
            injections(genericInfilename) = {};
        end
        injectionParFiles = injections(genericInfilename);
        injectionParFiles{sliceNum} = parfiles(i).name;
        injections(genericInfilename) = injectionParFiles;
    end

    
    %spread the joy of shifts
    injectionkeys = keys(injections);
    clear temp
    temp.jazz = 0;
    myholder = cell(length(injectionkeys),1);
    maxt = 0;
    for i=1:length(injectionkeys);
        parfiles = injections(injectionkeys{i});
        
        myholder{i} = cell(length(parfiles),1);
        for j=1:length(parfiles)
            if(~isempty(parfiles{j}))
                ParFileName = parfiles{j}; %#ok<NASGU>
                ReadPar;
                seriesNum = studyNum;
                if(~exist([outpath 'cinemri1.study' int2str(seriesNum) '.slice' int2str(sliceNum) '.mat']))
                    disp(['Could not find ' outpath 'cinemri1.study' int2str(seriesNum) '.slice' int2str(sliceNum) '.mat']);
                    continue;
                end
                disp(['Adding ' int2str(seriesNum) ' : ' int2str(sliceNum)]);
                load([outpath 'cinemri1.study' int2str(seriesNum) '.slice' int2str(sliceNum) '.mat']);

                [~,~,st] = size(cinemri1);
                maxt = max(maxt,st);
                myholder{i}{j}.cinemri1 = cinemri1;
                myholder{i}{j}.maxt = st;
                myholder{i}{j}.seriesNum = seriesNum;
                myholder{i}{j}.sliceNum = sliceNum;
                temp = strfind(injectionkeys{i},'/');
                temp = injectionkeys{i}((temp(end)+1):end);
                myholder{i}{j}.genericInjectionName = strrep(temp,'_','\_');
                myholder{i}{j}.blood = load([outpath 'blood_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)]);
                myholder{i}{j}.endo = load([outpath 'endo_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)]);
                myholder{i}{j}.epi = load([outpath 'epi_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)]);
                myholder{i}{j}.angle = load([outpath 'Roi_start_angle.study' int2str(seriesNum) '.slice' int2str(sliceNum)]);
            end
        end
    end
    figure(111);clf;
    a = size(myholder,1);
    b = 1;
    for i=1:a
        b = max(b,size(myholder{1},1));
    end
    for t=1:maxt
        for i=1:size(myholder,1)
            for j=1:size(myholder{i},1)
                if(~isempty(myholder{i}{j}))
                    h = subplot(b,a,(j-1)*a + i); imagesc(myholder{i}{j}.cinemri1(:,:,min(myholder{i}{j}.maxt,t))), colormap gray, hold on;
                    set(h,'XTick',[]);set(h,'YTick',[]);title(['Series ' num2str(myholder{i}{j}.seriesNum) ', Slice ' num2str(myholder{i}{j}.sliceNum)]);
                    plot(myholder{i}{j}.blood(:,1),myholder{i}{j}.blood(:,2));
                    plot(myholder{i}{j}.endo(:,1),myholder{i}{j}.endo(:,2));
                    plot(myholder{i}{j}.epi(:,1),myholder{i}{j}.epi(:,2));
                    if(j==1)
                        temp = {myholder{i}{j}.genericInjectionName,['Series ' num2str(myholder{i}{j}.seriesNum) ', Slice ' num2str(myholder{i}{j}.sliceNum)]};
                        if(mod(i,2) == 0)
                            temp = horzcat(temp{1},{''},temp{2});
                        end
                        title(temp);
                    end
                end
            end
        end
        pause();
    end