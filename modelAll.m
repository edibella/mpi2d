%function [cinemri1,time,shifts, endo,epi,blood,angle,curves,deltaSIcurves,AIF,pixelwiseUpslope,ktransValues,ve] = modelAll(pathToDicoms,writeFiles)
stagesToDo = [1 2];
%hard coded path to the dicom folder
pathToDicoms = '/v/raid1/bmatthew/MRIdata/Cardiac/Verio/P071610/ReconData/CV_Radial7Off_flexICE_fsa_newILOT32setupSTRESS4.2ML(41000)';
writeFiles = 0;

mypwd = pwd;
global waith
if(~exist('waith')),waith = 0;end
if(ishandle(waith))
    waith = waitbar(0,waith,'Good Morning Dave');
else
    waith = waitbar(0,'Good Morning Dave');
end
    

if(any(stagesToDo == 1))
    % set up variables
    cd(pathToDicoms);
    dicoms = dir('*.dcm');
    [sx sy] = size(dicomread(dicoms(1).name));
    st = length(dicoms);
    %% read in dicoms and upsample
    nonUpsampled = zeros(sx,sy,st);
    time = zeros(st,1);
    waitbar(0,waith,'Opening up Dicoms');
    header = dicominfo(dicoms(ceil(st/2)).name);
    figure(2), clf;
    figure(1),clf,subplot(2,2,1),cla, imagesc(dicomread(dicoms(1).name)), colormap gray, title(strrep(header.SeriesDescription,'_','\_')),pause(.1);
    for i=1:length(dicoms)
        nonUpsampled(:,:,i) = dicomread(dicoms(i).name);
        header = dicominfo(dicoms(i).name);
        timestr = header.AcquisitionTime;
        time(i) = str2num(timestr(1:2))*60*60 + str2num(timestr(3:4))*60 + str2num(timestr(5:end));
        seriesNum = header.SeriesNumber;
        waitbar(i/st,waith,num2str(time(i)));
    end
    cd(mypwd);   %go back
    [time IX] = sort(time);
    time = time-time(1);
    nonUpsampled = nonUpsampled(:,:,IX);
    %find the slice number
    temp = pathToDicoms(1:strfind(pathToDicoms,'('));
    base = pathToDicoms(strfind(pathToDicoms,'/'):end);
    results = dir([temp '*']);
    sliceNum = 0;
    for i=1:length(results)
        sliceNum = sliceNum + 1;
        if(strcmp(results,base))
            break
        end
    end
    %% upsample
    temp = interp2(nonUpsampled(:,:,1),2,'cubic');
    [sx sy] = size(temp);
    fullcinemri1 = zeros(sx,sy,st);
    waitbar(0,waith,'Upsampling');
    for t=1:st
        waitbar(t/st,waith);
        fullcinemri1(:,:,t) = interp2(nonUpsampled(:,:,t),2,'cubic');
    end
    %% pick the LVRV
    [RV,LV] = FindLVRV(fullcinemri1,0);
    rangex = (LV(1) - 75):(LV(1) + 75);
    rangey = (LV(2) - 75):(LV(2) + 75);
    RV(1) = max(RV(1),min(rangex)+1);
    RV(2) = max(RV(2),min(rangey)+1);
    temp = zeros(size(fullcinemri1,1),size(fullcinemri1,2),3);
    temp2 = fullcinemri1(:,:,22);
    temp2 = temp2-min(temp2(:));
    temp2 = temp2/max(temp2(:));
    temp(:,:,1) = temp2;
    temp(:,:,2) = temp2;
    temp(:,:,3) = temp2;
    figure(1),subplot(2,2,1), imagesc(temp), colormap gray, hold on;
    title([strrep(header.SeriesDescription,'_','\_') ' and cropping']);
    plot(LV(2), LV(1),'ro',RV(2), RV(1),'o');
    text(LV(2)+10, LV(1)-15,'LV','Color',[1 1 1]);
    text(RV(2)+10, RV(1)-15,'RV','Color',[1 1 1]);
    plot([min(rangey) max(rangey) max(rangey) min(rangey) min(rangey)],[min(rangex) min(rangex) max(rangex) max(rangex) min(rangex)]);
    RV = RV - [min(rangex) min(rangey)];
    LV = LV - [min(rangex) min(rangey)];
    if(writeFiles)
        %create the par file
        ParFileName = ['series' num2str(seriesNum) '.slice' num2str(sliceNum) '.par'];
        copyfile('Template.par',ParFileName);
        clear temp;  
        temp.rangex = [num2str(min(rangex)) ':' num2str(max(rangex))];
        temp.rangey = [num2str(min(rangey)) ':' num2str(max(rangey))];
        temp.infile = [pathToDicoms '/'];
        temp.studyNum = seriesNum;
        temp.sliceNum = sliceNum;
        temp.seriesNumAIF = seriesNum;
        temp.sliceNumAIF = sliceNum;
        temp.ranget = ['2:' num2str(st)];
        temp.lastFrame = st;
        temp.timeStampFile = ['timeStampSer' num2str(seriesNum) '.mat'];
        updateParFile(ParFileName,temp);
        
        %create the timestamp file
        timeStamp = time;
        save(temp.timeStampFile,'timeStamp');
        
        save(['Output/RVLVPos.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.mat'],'RV','LV');
    end


    %% modelBased Registration
    
    [shifts,parameters,delays,pixelwisevariance,cinemri1] = modelRegistration(fullcinemri1,rangex, rangey,0,2);
    subplot(2,2,2), cla,hold on
    plot(shifts(:,1));
    plot(shifts(:,2)+5);
    title('Shifts(Top:Y, Bottom:X)');
    subplot(2,2,3), imagesc(delays), colorbar, title('pixelwise Delays(time @ onset)');
    mycine = zeros(length(rangex),length(rangey),st);
    for t=1:length(shifts)
        temp = circshift(fullcinemri1(:,:,t),shifts(t,:));
        mycine(:,:,t) = temp(rangex,rangey);
    end
    cinemri1 = mycine;
    
    if(writeFiles)
        save(['Output/cinemri1.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.mat'],'cinemri1');
        fid = fopen(['Output/shiftsMAN.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.txt']);
        for t=1:length(shifts)
            fprintf(fid,'%f	%f\n',shifts(t,:));
        end
        fclose(fid);
        
        save(['Output/Variance.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.mat'],'pixelwisevariance');
        save(['Output/Delays.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.mat'],'delays');
    end
end

%% model based Segmentation
if(any(stagesToDo == 2))
    %variables needed parameters, delays, cinemri1, RV, LV, 
    
    [sx sy st] = size(cinemri1);
    clear smoothedMovie
    mymean = mean(cinemri1(:,:,1:5),3);
    for t=1:st
        cinemri1(:,:,t) = cinemri1(:,:,t) - mymean;
    end

    %% make a temporally smooth curve for each pixel
    side = 5;
    st2 = floor(st);
    myvals = zeros(sx*sy,length((1+side):st2));
    padded = cat(3,repmat(mean(cinemri1(:,:,1:side),3),[1 1 side]), cinemri1);
    h = diff(fspecial('gaussian',[st 1],10));
    temp = zeros(sx,sy,length((1+side):st2));
    for i=1:sx
        waitbar(i/sx,waith,'Creating Temporally smooth movie');
        for j=1:sy
            curve = squeeze(padded(i,j,:));
            
            curve = filter2(h,curve);
            myvals((j-1)*sx+i,:) = curve((1+side):st2);
            smoothedMovie(i,j,:) = curve((1+side):st2);
            %convert to deltaSI
        end
    end

    %initialize 4 clusters
    p1 = squeeze(padded(RV(1),RV(2),(1+side):st2))';
    p2 = squeeze(padded(LV(1),LV(2),(1+side):st2))';
    temp = round(mean(vertcat(RV,LV),1));
    p3 = squeeze(padded(temp(1),temp(2),(1+side):st2))';
    p4 = squeeze(padded(1,1,(1+side):st2))';

    %run kmeans and create a master mask
    %IDX = kmeans(myvals,4,'start',vertcat(p1,p2,p3,p4),'emptyaction','singleton');
    waitbar(0,waith,'K-means grouping of Temporally smooth movie');
    IDX = kmeans(myvals,4,'emptyaction','singleton');
    master = reshape(IDX,[sx sy]);
    figure(2)
    subplot(2,2,1),imagesc(master), colormap(lines(6));h = colorbar; title('K-means grouping');
    cbfreeze(h);
    LVMask = imfill(master== master(LV(1), LV(2)),LV);
    RVMask = imfill(master== master(RV(1), RV(2)),RV);
    LVperimiter = bwperim(LVMask);
    oldtemporalMaster = master;
    originalDelays = delays;
    mydelays = delays;
    
    b1 = parameters(:,:,1);%Efficient Method for Calculating Kinetic Parameters Using
    b2 = parameters(:,:,2);%T1-Weighted Dynamic Contrast-Enhanced Magnetic
    b3 = parameters(:,:,3);%Resonance Imaging

    
    
    %find myocardium by finding the space between LV and RV
    [rl,cl] = find(LVMask);
    meanLV = mean([rl cl]);
    [rr,cr] = find(RVMask);
    meanRV = mean([rr cr]);
    
    
    %find myocardium between LV and RV
    [thetar,radiusr] = cart2pol(rr-meanLV(1),cr-meanLV(2));
    [mintheta,minthetai] = min(thetar);
    [topRV(1), topRV(2)] = pol2cart(thetar(minthetai), radiusr(minthetai));
    topRV = topRV + meanLV;
    figure(1), subplot(2,2,1), 
    plot(topRV(2) + min(rangey), topRV(1) + min(rangex),'o');
    plot([meanLV(2) topRV(2)] + min(rangey),[meanLV(1) topRV(1)] + min(rangex),'r');
    %plot(topRV(2),topRV(1), 'o');
    [maxtheta,maxthetai] = max(thetar);
    [bottomRV(1), bottomRV(2)] = pol2cart(thetar(maxthetai), radiusr(maxthetai));
    bottomRV = bottomRV + meanLV;
    %plot(bottomRV(2),bottomRV(1), 'o');
    
    
    [thetal,radiusl] = cart2pol(rl-meanRV(1),cl-meanRV(2));
    [maxtheta,maxthetai] = max(thetal);
    [topLV(1), topLV(2)] = pol2cart(thetal(maxthetai), radiusl(maxthetai));
    topLV = topLV + meanRV;
    %plot(topLV(2),topLV(1), 'o');
    [mintheta,minthetai] = min(thetal);
    [bottomLV(1), bottomLV(2)] = pol2cart(thetal(minthetai), radiusl(minthetai));
    bottomLV = bottomLV + meanRV;
    %plot(bottomLV(2),bottomLV(1), 'o');
    
    prelimMyo = (roipoly(master,[topRV(2) meanRV(2) bottomRV(2) bottomLV(2) meanLV(2) topLV(2)],[topRV(1) meanRV(1) bottomRV(1) bottomLV(1) meanLV(1) topLV(1)]) - LVMask - RVMask)>0; 
    
    
    [r,c] = find(prelimMyo>0);
    [theta,radius] = cart2pol(r-meanLV(1),c-meanLV(2));
    meanRadius = mean(radius);
    fatness = std(radius);
    prelimCircle = vertcat(meanRadius*cos(0:.01:(2*pi))+meanLV(1), meanRadius*sin(0:.01:(2*pi))+meanLV(2));
    temp = zeros(sx,sy);
    for i=1:length(prelimCircle)
        if(~isnan(prelimCircle(1,i)))
            temp(ceil(prelimCircle(1,i)),ceil(prelimCircle(2,i))) = 1;
        end
    end
    minimumConnectivity = bwmorph(temp,'dilate',1);
    mygaussian = fspecial('gaussian',[sx,sy],4*fatness);
    mygaussian = mygaussian/mygaussian(round(sx/2),round(sy/2));
    mygaussian = mygaussian .*mygaussian.*mygaussian.*mygaussian;
    prelimMyo = filter2(mygaussian,temp);
    prelimMyoWeight = prelimMyo / mean(prelimMyo(temp>0));
    prelimMyo = ((prelimMyoWeight > .5) - LVMask - RVMask)>0;
    
    
    prelimMyo1 = (bwmorph(LVMask,'dilate',15) - LVMask - RVMask) > 0;
    if(sum(sum(prelimMyo1.*prelimMyo))/sum(sum(prelimMyo1)) < .8)
        prelimMyo = prelimMyo1;
    else
        prelimMyo = (prelimMyo + prelimMyo1) > 0;
    end
    
    
    %figure(5), imagesc(prelimMyo), colorbar;
    
    
    
    myslopefinder = -diff(fspecial('gaussian',[(st+1),1],7));
    pixelwiseUpslope = zeros(sx,sy);
    for x=1:sx
        for y=1:sy
            if(prelimMyo(x,y) > 0)
                pixelwiseUpslope(x,y) = max(filter2(myslopefinder,squeeze(mycine(x,y,:))));
            end
        end
    end
    
    
    %figure(2)
    %subplot(2,2,1),imagesc(master), colorbar, title('Master regions');
    %subplot(2,2,2),imagesc(mydelays),colorbar, title('Delays');
    %subplot(2,2,3), imagesc(ktrans), colorbar, title('ktrans Master');
    
    
    upslopemean = median(pixelwiseUpslope(prelimMyo>0));
    upslopestd = std(pixelwiseUpslope(prelimMyo>0));
    
    variancemean = median(pixelwisevariance(prelimMyo>0));
    variancestd = std(pixelwisevariance(prelimMyo>0));
    
    delaymean = median(mydelays(prelimMyo>0));
    delaystd = std(mydelays(prelimMyo>0));
    
    b1mean = median(b1(prelimMyo>0));
    b1std = std(b1(prelimMyo>0));
    b2mean = median(b2(prelimMyo>0));
    b2std = std(b2(prelimMyo>0));
    b3mean = median(b3(prelimMyo>0));
    b3std = std(b3(prelimMyo>0));
    
    
    %temporalKmeans' contribution
    masterMask = master == median(master(prelimMyo>0));
    delayMask = (mydelays < (delaymean + 1.5*delaystd)) .* (mydelays > (delaymean - 1.5*delaystd));
    b1Mask = (b1 < (b1mean + 1.5*b1std)) .* (b1 > (b1mean - 1.5*b1std));
    b2Mask = (b2 < (b2mean + 1.5*b2std)) .* (b2 > (b2mean - 1.5*b2std));
    b3Mask = (b3 < (b3mean + 1.5*b3std)) .* (b3 > (b3mean - 1.5*b3std));
    varianceMask = (b3 < (b3mean + 1.5*b3std)) .* (b3 > (b3mean - 1.5*b3std));
    upslopeMask = (pixelwiseUpslope < (upslopemean + 1.5*upslopestd)) .* (pixelwiseUpslope > (upslopemean - 1.5*upslopestd));
    
    delayWeigth = 1-1./(1+(mydelays/delaymean).^6);
    b1Weigth = exp(-(b1 -b1mean).^2 /(2*(b1std)^2));
    b2Weigth = exp(-(b2 -b2mean).^2 /(2*(b2std)^2));
    b3Weigth = exp(-(b3 -b2mean).^2 /(2*(b3std)^2));
    varianceWeight = exp(-(pixelwisevariance -variancemean).^2 /(2*(variancestd)^2));
    upslopeWeight = exp(-(pixelwiseUpslope -upslopemean).^2 /(2*(upslopestd)^2));
    
    %figure(4)
    %subplot(2,2,1), imagesc(delayMask), title('delay mask');
    %subplot(2,2,2), imagesc(ktransMask), title('ktrans mask');
    %subplot(2,2,3), imagesc(masterMask), title('master mask');
    if(mean(delayMask(prelimMyo>0)) > .5)
        prelimMyo = prelimMyo .* delayMask;
        prelimMyoWeight = prelimMyoWeight .* delayWeigth;
    end
    if(mean(b1Mask(prelimMyo>0)) > .5)
        prelimMyo = prelimMyo .* b1Mask;
        prelimMyoWeight = prelimMyoWeight .* b1Weigth;
    end
    if(mean(b2Mask(prelimMyo>0)) > .5)
        prelimMyo = prelimMyo .* b2Mask;
        prelimMyoWeight = prelimMyoWeight .* b2Weigth;
    end
    if(mean(b3Mask(prelimMyo>0)) > .5)
        prelimMyo = prelimMyo .* b3Mask;
        prelimMyoWeight = prelimMyoWeight .* b3Weigth;
    end
    if(mean(masterMask(prelimMyo>0)) > .5)
        prelimMyo = prelimMyo .* masterMask;
        prelimMyoWeight = prelimMyoWeight .* masterMask;
    end
    if(mean(varianceMask(prelimMyo>0)) > .5)
        prelimMyo = prelimMyo .* varianceMask;
        prelimMyoWeight = prelimMyoWeight .* varianceWeight;
    end
    if(mean(upslopeMask(prelimMyo>0)) > .5)
        prelimMyo = prelimMyo .* upslopeMask;
        prelimMyoWeight = prelimMyoWeight .* upslopeWeight;
    end
    mygaussian = fspecial('gaussian',[sx,sy],2);
    %myo = filter2(mygaussian,prelimMyo)>.5;
    prelimMyoWeight = prelimMyoWeight / median(prelimMyoWeight(prelimMyoWeight>0));
    myo = prelimMyoWeight>.5;
    figure(2), subplot(2,2,2), imagesc(prelimMyoWeight), colorbar, title('Weighted myocardium, threshold@ .01');
    
    %ensure it's connected
    %[myo,minimumConnectivity] = maintainMinumumConnectivity(myo,minimumConnectivity);
    
    myinner = imfill(~myo,LV);
    myinner = filter2(fspecial('gaussian',[sx,sy],5),myinner)>.5;
    myo = (myo - myinner)>0;
    
    %smooth the myocardium
    myo = filter2(fspecial('gaussian',[sx,sy],3),myo)>.5;
    
    
    
%     mygaussian = fspecial('gaussian',[sx sy],40);
%     weightingTerm = filter2(mygaussian,LVperimiter);
%     weightingTerm = weightingTerm / median(weightingTerm(LVperimiter>0));
%     weightingTerm = prelimMyoWeight;
%     mastersmyo = master== median(master(prelimMyo>0));
%     ktransmyo = ktransMaster == median(ktransMaster(mastersmyo>0));
%     myvariances = pixelwisevariance(ktransmyo>0);
%     medianVariance = median(myvariances);
%     stdVariance = std(myvariances);
%     varianceMyo = (pixelwisevariance < medianVariance + 2*stdVariance).*ktransmyo;
%     
%     mygaussian = fspecial('gaussian',[sx sy],10);
%     myo = (filter2(mygaussian,mydelays)+filter2(2*mygaussian,ktransmyo)+filter2(mygaussian,mastersmyo)+filter2(mygaussian,varianceMyo)).*weightingTerm>1.2;
%     myo = myo - LVMask;
%     subplot(2,2,4), imagesc(myo), colorbar, title('final');
    
    
    
    
    figure(2), subplot(2,2,3);
    myrgb = zeros(sx,sy,3);
    temp = cinemri1(:,:,22);
    temp = temp - min(temp(:));
    temp = temp/max(temp(:));
    myrgb(:,:,2) = temp;
    myrgb(:,:,3) = temp;
    temp(myo > 0) = min(1,temp(myo > 0)*2);
    myrgb(:,:,1) = temp;
    imagesc(myrgb);
    
    %now make contours that extracts out the myocardium
    [r,c] = find(myo>0);
    [theta,radius] = cart2pol(r-meanLV(1),c- meanLV(2));
    meanR = mean(radius);
    connectedMyo = myo;
    for theta = 0:.01:(2*pi)
        x = round(meanR*cos(theta) + meanLV(1));
        y = round(meanR*sin(theta) + meanLV(2));
        connectedMyo((x-1):(x+1),(y-1):(y+1)) = 1;
    end
    %endocardium
    inner = imfill(~connectedMyo,round(meanLV));
    [r,c] = find(inner>0);
    [minr, minri] = min(r);
    boundary = bwtraceboundary(inner,[r(minri) c(minri)],'N');
    [theta,r] = cart2pol(boundary(:,1) -meanLV(1), boundary(:,2) -meanLV(2));
    r = r * 1.01;
    [x,y] = pol2cart(theta,r);
    x = x + meanLV(1); y = y + meanLV(2);
    endo = horzcat(x,y);
    
    %epicardium
    endoMask = connectedMyo + inner;
    [r,c] = find(endoMask>0);
    [minr, minri] = min(r);
    boundary = bwtraceboundary(endoMask,[r(minri) c(minri)],'N');
    [theta,r] = cart2pol(boundary(:,1) -meanLV(1), boundary(:,2) -meanLV(2));
    r = r * 1.01;
    [x,y] = pol2cart(theta,r);
    x = x + meanLV(1); y = y + meanLV(2);
    epi = horzcat(x,y);
    
    %find the contour for the blood pool that gives the highest results
%     [r,c] = find(LVMask);
%     bestx = r(1); besty = c(1);
%     highestPeak = max(mean(mean(cinemri1((bestx-3):(bestx+3),(besty-3):(besty+3),:),1),2));
%     for i=1:length(r)
%         if(max(mean(mean(cinemri1((r(i)-3):(r(i)+3),(c(i)-3):(c(i)+3),:),1),2)) > highestPeak)
%             bestx = r(i); besty = c(i);
%             highestPeak = max(mean(mean(cinemri1((r(i)-3):(r(i)+3),(c(i)-3):(c(i)+3),:),1),2));
%         end
%     end
%     LVminiMask = zeros(sx,sy);
%     LVminiMask((bestx-3):(bestx+3),(besty-3):(besty+3)) = 1;
%     [r,c] = find(LVminiMask>0);
%     [minr, minri] = min(r);
%     boundary = bwtraceboundary(LVminiMask,[r(minri) c(minri)],'N');
%     [theta,r] = cart2pol(boundary(:,1) -meanLV(1), boundary(:,2) -meanLV(2));
%     r = r * 1.01;
%     [x,y] = pol2cart(theta,r);
%     x = x + meanLV(1); y = y + meanLV(2);
%     blood = horzcat(x,y);
    
    
    [r,c] = find(LVMask>0);
    [minr, minri] = min(r);
    boundary = bwtraceboundary(LVMask,[r(minri) c(minri)],'N');
    [theta,r] = cart2pol(boundary(:,1) -meanLV(1), boundary(:,2) -meanLV(2));
    r = r * 1.01;
    [x,y] = pol2cart(theta,r);
    x = x + meanLV(1); y = y + meanLV(2);
    blood = horzcat(x,y);
    
    %find the angle to the LV RV insertion point
    [angle,radius] = cart2pol(topRV(1) - meanLV(1),topRV(2) - meanLV(2));
    
    %use the contours to recreate a mask as that's what the user would be
    %doing to recreate our results
    outtermask = roipoly(cinemri1(:,:,22),epi(:,2),epi(:,1));
    innermask = roipoly(cinemri1(:,:,22),endo(:,2),endo(:,1));
    newMyo = (outtermask - innermask)>0;
    newLV = roipoly(cinemri1(:,:,22),blood(:,2),blood(:,1));
    [r,c] = find(newLV>0);
    newMeanLV = [mean(r) mean(c)];
    
    %use the myocardial mask to find the curves then the deltaSI curves
    [r,c] = find(outtermask>0);
    middle = [mean(r) mean(c)];
    [r,c] = find(newMyo>0);
    [theta,radius] = cart2pol(r-middle(1), c - middle(2));
    theta = mod(theta + angle,2*pi);
    
    mycolors = reshape(lines(6),1,6,3);
    figure(2), subplot(2,2,3);
    myrgb = zeros(sx,sy,3);
    temp = cinemri1(:,:,22);
    temp = temp - min(temp(:));
    temp = temp/max(temp(:));
    myrgb(:,:,2) = temp;
    myrgb(:,:,3) = temp;
    myrgb(:,:,1) = temp;
    myregionalMask = zeros(sx,sy);
    labeled = floor(theta/(2*pi)*6);
    for i=1:length(r)
        myregionalMask(r(i),c(i)) = labeled(i)+1;
        myrgb(r(i),c(i),:) = myrgb(r(i),c(i),:).*mycolors(1,labeled(i)+1,:);
    end
    imagesc(myrgb);
    deltaSIconverter = median(cinemri1(:,:,1:7),3);
    mycine = cinemri1;
    for t=1:st
        mycine(:,:,t) = mycine(:,:,t) - deltaSIconverter;
    end
    deltaSIcurves = zeros(st,6);
    mediandeltaSIcurves = zeros(st,6);
    curves = zeros(st,6);
    for region=1:6
        regionalMask = myregionalMask == region;
        for t=1:st
            temp = mycine(:,:,t);
            temp2 = cinemri1(:,:,t);
            curves(t,region) = mean(temp2(regionalMask>0));
            deltaSIcurves(t,region) = mean(temp(regionalMask>0));
            mediandeltaSIcurves(t,region) = median(temp(regionalMask>0));
        end
    end
    AIF = zeros(st,1);
    medianAIF = zeros(st,1);
    fifthHighestInAIF = zeros(st,1);
    
    for t=1:st
        temp = mycine(:,:,t);
        values = temp(newLV>0);
        values = sort(values);
        AIF(t) = mean(temp(newLV>0));
        medianAIF(t) = median(temp(newLV>0));
        fifthHighestInAIF(t) = values(end-5);
    end
    AIF = AIF - mean(AIF(1:5));
    medianAIF = medianAIF - mean(medianAIF(1:5));
    fifthHighestInAIF = fifthHighestInAIF - mean(fifthHighestInAIF(1:5));
    %do pixelwise ratio of pixel's max upslope / AIF max upslope
    myslopefinder = -diff(fspecial('gaussian',[(st+1),1],7));
    AIFupslope = max(filter2(myslopefinder,fifthHighestInAIF));
    pixelwiseUpslope = zeros(sx,sy);
    for x=1:sx
        for y=1:sy
            if(myregionalMask(x,y) > 0)
                pixelwiseUpslope(x,y) = max(filter2(myslopefinder,squeeze(mycine(x,y,:))))/AIFupslope;
            end
        end
    end
    smallest = min(pixelwiseUpslope(pixelwiseUpslope > 0));
    temp = pixelwiseUpslope;
    temp(temp == 0) = .8*smallest;
    figure(2), subplot(2,2,2), imagesc(temp), title('pixelwise ratio of max upslope over AIF upslope'), colorbar;
    %figure(4),clf, hold on, title('average deltaSI Curves');
    figure(2), subplot(2,2,3), cla
    figure(2), subplot(2,2,4),cla, hold on, title('median deltaSI Curves');
    runningMask = zeros(sx,sy);
    for region=1:6
        regionalMask = myregionalMask == region;
        runningMask = runningMask + regionalMask;
        %figure(4)
        %plot(deltaSIcurves(:,region),'Color',mycolors(1,region,:));
        figure(2), subplot(2,2,4),
        plot(mediandeltaSIcurves(:,region),'Color',mycolors(1,region,:));
        figure(2), subplot(2,2,3),
        temp = myrgb;
        temp(:,:,1) = temp(:,:,1) .* runningMask;
        temp(:,:,2) = temp(:,:,2) .* runningMask;
        temp(:,:,3) = temp(:,:,3) .* runningMask;
        imagesc(temp);
        pause(.5);
    end
    figure(2), subplot(2,2,3);
    imagesc(myrgb);
    %figure(4);
    %plot(AIF,'k');
    figure(2), subplot(2,2,4),
    plot(fifthHighestInAIF,'k');
    AIF = fifthHighestInAIF;
    
    if(writeFiles)
        %save the contours, curves/deltaSIcurves + AIF, pixelwiseUpslope
        fid = fopen(['Output/blood_polyCoords.study' num2str(seriesNum) '.slice' num2str(sliceNum)]);
        for t=1:length(blood)
            fprintf(fid,'%f	%f\n',blood(t,:));
        end
        fclose(fid);
        fid = fopen(['Output/endo_polyCoords.study' num2str(seriesNum) '.slice' num2str(sliceNum)]);
        for t=1:length(endo)
            fprintf(fid,'%f	%f\n',endo(t,:));
        end
        fclose(fid);
        fid = fopen(['Output/epi_polyCoords.study' num2str(seriesNum) '.slice' num2str(sliceNum)]);
        for t=1:length(epi)
            fprintf(fid,'%f	%f\n',epi(t,:));
        end
        fclose(fid);
        fid = fopen(['Output/Roi_start_angle.study' num2str(seriesNum) '.slice' num2str(sliceNum)]);
        fprintf(fid,'%f\n',angle);
        fclose(fid);
        
        temp = curves;
        curves = vertcat(AIF',curves');
        save(['Output/curves.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'curves');
        curves = temp;
        
        deltaSIcurves = vertcat(AIF',mediandeltaSIcurves');
        save(['Output/deltaSIcurves.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'deltaSIcurves');
        upslope = pixelwiseUpslope;
        save(['Output/upslope.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'upslope');
        
        save(['Output/upslope.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'upslope');
    end
end

if(any(stagesToDo == 3))
    %variables needed: deltaSIcurves, seriesNum, sliceNum, AIF
    %save deltaSIcurves to disk and run mpi2d stg=4.1
    cd(mypwd);
    template = ['.study' num2str(seriesNum) '.slice' num2str(sliceNum)];
    deltaSIcurves = vertcat(AIF',mediandeltaSIcurves');
    
    save(['deltaSIcurves' template '.mat'],'deltaSIcurves');
    mpi2d stg=4.1
    flowFiles = dir(['Output/flowvalues' template '.6.1_fixedDelay0.AIF_' num2str(seriesNum) '_' num2str(sliceNum) '_1.txt.full.txt']);
    if(length(flowFiles) == 0)
        x = 0;
    else
        clear ktrans ve
        ktrans = zeros(6,1);
        ve = zeros(6,1);
        fid = fopen(flowFiles);
        tline = fgetl(fid);
        while(ischar(tline))
            if(~isempty(strfind(tline,'Ktrans')))
                A = sscanf(tline,'Ktrans  %d     %f');
                ktrans(A(1)) = A(2);
            end
            if(~isempty(strfind(tline,'ve')))
                A = sscanf(tline,'ve    %d      %f');
                ve(A(1)) = A(2);
            end
            tline = fgetl(fid);
        end
        fclose(fid);
    end
    if(~writeFiles)
        delete(['Output/' flowFiles(1).name]);
    end
end
