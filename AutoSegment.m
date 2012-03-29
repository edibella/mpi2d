function [blood ,epi,endo,angle] = AutoSegment(seriesNum,sliceNum,myothreshold,UseExistingAsPreliminary)
    if(~exist('UseExistingAsPreliminary','var'))
        UseExistingAsPreliminary = 0;
    end
    blood = [];
    epi = [];
    endo = [];
    angle = [];
    infilename=strcat('Output/','cinemri1.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.mat');
    load(infilename);
    %construct pixelwise variance, delays, k-trans images
    mask = ones(size(cinemri1,1),size(cinemri1,2));
    if(exist(['Output/RVLVpos.study' int2str(seriesNum) '.slice' int2str(sliceNum) '.mat']))
        load(['Output/RVLVpos.study' int2str(seriesNum) '.slice' int2str(sliceNum) '.mat']);
    else
        [RV,LV] = FindLVRV(cinemri1);
    end
    %RVcurve = squeeze(mean(mean(cinemri1((RV(1)-2):(RV(1)+2),(RV(2)-2):(RV(2)+2),:),1),2));
    LVcurve = squeeze(mean(mean(cinemri1((LV(1)-2):(LV(1)+2),(LV(2)-2):(LV(2)+2),:),1),2));
    [~,LVonset] = max(LVcurve);
    if(LVonset+16 > size(cinemri1,3))
        figure(1), clf, plot(LVcurve);
        disp('If this looks like a proper AIF then type dbcontinue.  If it doesn''t please find an LV point that does and type dbcontinue.  If this dataset is unsatisfactory then type dbquit');
        %keyboard;   %bmatthew we should be asking the user if it's an OK LV point, but this algorithm should do fine without a proper LV curve
        LVcurve = squeeze(mean(mean(cinemri1((LV(1)-2):(LV(1)+2),(LV(2)-2):(LV(2)+2),:),1),2));
        [~,LVonset] = max(LVcurve);
    end
    TimeDelayWindow = max(1,round(.8*LVonset)):min(round(LVonset+16),size(cinemri1,3));

    global waith %#ok<TLEV>
    if(~exist('waith','var')),waith = 0;end
    if(ishandle(waith))
        waith = waitbar(0,waith,'Good Morning Dave');
    else
        waith = waitbar(0,'Good Morning Dave');
    end
    waitbar(0,waith,'Extracting model params');
    [parameters,delays,smoothed] = model(squeeze(LVcurve),TimeDelayWindow,cinemri1,mask); %#ok<NASGU>
  
    figure(seriesNum*10 + 2*sliceNum)
    subplot(2,2,3), cla;
    imagesc(delays), colormap jet, colorbar, title('Delays');
    
    
    
    [sx sy st] = size(cinemri1);
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
    for i=1:sx
        waitbar(i/sx,waith,'Creating Temporally smooth movie');
        for j=1:sy
            curve = squeeze(padded(i,j,:));
            curve = filter2(h,curve);
            myvals((j-1)*sx+i,:) = curve((1+side):st2);
            mydeltaSImovie(i,j,:) = curve((1+side):st2);
        end
    end

    %initialize 4 clusters
    p1 = squeeze(padded(RV(1),RV(2),(1+side):st2))';
    p2 = squeeze(padded(LV(1),LV(2),(1+side):st2))';
    temp = round(mean(vertcat(RV,LV),1));
    p3 = squeeze(padded(temp(1),temp(2),(1+side):st2))';
    p4 = squeeze(padded(1,1,(1+side):st2))';

    %run kmeans and create a master mask
    waitbar(0,waith,'K-means grouping of Temporally smooth movie');
    IDX = kmeans(myvals,4,'start',vertcat(p1,p2,p3,p4),'emptyaction','singleton');
    %IDX = kmeans(myvals,4,'emptyaction','singleton');
    master = reshape(IDX,[sx sy]);
    %figure(seriesNum*10 + 2*sliceNum+1), subplot(2,2,1)
    %imagesc(master), colormap(lines(numAzimuthalRegions));h = colorbar; title('K-means grouping');cbfreeze(h);
    close(waith);
    LVMask = imfill(master== master(LV(1), LV(2)),LV);
    RVMask = imfill(master== master(RV(1), RV(2)),RV);
    mydelays = delays;
    
    b1 = parameters(:,:,1);%Efficient Method for Calculating Kinetic Parameters Using
    b2 = parameters(:,:,2);%T1-Weighted Dynamic Contrast-Enhanced Magnetic
    b3 = parameters(:,:,3);%Resonance Imaging

    
    %find myocardium by finding the space between LV and RV
    [rl,cl] = find(LVMask);
    meanLV = mean([rl cl]);
    [rr,cr] = find(RVMask);
    meanRV = mean([rr cr]);
    if(~UseExistingAsPreliminary)


        %find myocardium between LV and RV
        [thetar,radiusr] = cart2pol(rr-meanLV(1),cr-meanLV(2));
        [~,minthetai] = min(thetar);
        [topRV(1), topRV(2)] = pol2cart(thetar(minthetai), radiusr(minthetai));
        topRV = topRV + meanLV;
        %temp = cinemri1(:,:,22);
        %temp = temp - min(temp(:));
        %temp = temp/max(temp(:));
        %myrgb(:,:,1) = temp;
        %myrgb(:,:,2) = temp;
        %myrgb(:,:,3) = temp;
        %figure(seriesNum*10 + 2*sliceNum),subplot(2,2,1), cla; imagesc(myrgb), hold on,
        %plot(meanLV(2), meanLV(1),'o');
        %plot(meanRV(2), meanRV(1),'ro');
        %plot([meanLV(2) topRV(2)],[meanLV(1) topRV(1)],'r');
        %plot(topRV(2),topRV(1), 'o');
        [~,maxthetai] = max(thetar);
        [bottomRV(1), bottomRV(2)] = pol2cart(thetar(maxthetai), radiusr(maxthetai));
        bottomRV = bottomRV + meanLV;
        %plot(bottomRV(2),bottomRV(1), 'o');


        [thetal,radiusl] = cart2pol(rl-meanRV(1),cl-meanRV(2));
        [~,maxthetai] = max(thetal);
        [topLV(1), topLV(2)] = pol2cart(thetal(maxthetai), radiusl(maxthetai));
        topLV = topLV + meanRV;
        %plot(topLV(2),topLV(1), 'o');
        [~,minthetai] = min(thetal);
        [bottomLV(1), bottomLV(2)] = pol2cart(thetal(minthetai), radiusl(minthetai));
        bottomLV = bottomLV + meanRV;
        %plot(bottomLV(2),bottomLV(1), 'o');

        prelimMyo = (roipoly(master,[topRV(2) meanRV(2) bottomRV(2) bottomLV(2) meanLV(2) topLV(2)],[topRV(1) meanRV(1) bottomRV(1) bottomLV(1) meanLV(1) topLV(1)]) - LVMask)>0; 

        if(sum(prelimMyo(:)) > 100)

            [r,c] = find(prelimMyo>0);
            [~,radius] = cart2pol(r-meanLV(1),c-meanLV(2));
            meanRadius = mean(radius);
            fatness = std(radius);
            prelimCircle = vertcat(meanRadius*cos(0:.01:(2*pi))+meanLV(1), meanRadius*sin(0:.01:(2*pi))+meanLV(2));
            %saturate circle
            prelimCircle(prelimCircle < 0) = 1;
            prelimCircle(1,prelimCircle(1,:) > sx) = sx;
            prelimCircle(2,prelimCircle(2,:) > sy) = sy;
            temp = zeros(sx,sy);
            for i=1:length(prelimCircle)
                if(~isnan(prelimCircle(1,i)))
                    temp(ceil(prelimCircle(1,i)),ceil(prelimCircle(2,i))) = 1;
                end
            end
            mygaussian = fspecial('gaussian',[sx,sy],4*fatness);
            mygaussian = mygaussian/mygaussian(round(sx/2),round(sy/2));
            mygaussian = mygaussian .*mygaussian.*mygaussian.*mygaussian;
            prelimMyo = filter2(mygaussian,temp);
            prelimMyoWeight = prelimMyo / mean(prelimMyo(temp>0));
            prelimMyo = ((prelimMyoWeight > .5) - bwmorph(LVMask,'dilate',4))>0;
            if(median(master(prelimMyo>0)) ~= master(1,1))
                prelimMyo = ((prelimMyoWeight > .5) - bwmorph(LVMask,'dilate',4) - (master==master(1,1)))>0;
            else
                prelimMyo = ((prelimMyoWeight > .5) - bwmorph(LVMask,'dilate',4))>0;
            end

        else
            prelimMyo = (bwmorph(LVMask,'dilate',sx/8) - LVMask-(master==master(1,1))) > 0;
            prelimMyoWeight = prelimMyo;
        end
    else
        endofile = ['Output/endo_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)];
        epifile = ['Output/epi_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)];
        endo = load(endofile);
        epi = load(epifile);
        inner = roipoly(cinemri1(:,:,22), endo(:,1),endo(:,2));
        outer = roipoly(cinemri1(:,:,22), epi(:,1),epi(:,2));
        prelimMyo = (outer - inner) > 0;
        prelimMyoWeight = bwmorph(prelimMyo,'dilate',5);
    end
    
    
    
    
    
    
    %like k-means, but I'm getting the distance to the average
    %prelimMyotissue curve
    tst = size(mydeltaSImovie,3);
    meanTissue = zeros(tst,1);
    for t=1:tst
        temp = mydeltaSImovie(:,:,t);
        meanTissue(t) = mean(temp(prelimMyo));
    end
    mydistimage = zeros(sx,sy);
    for x=1:sx
        for y=1:sy
            mydistimage(x,y) = norm(squeeze(mydeltaSImovie(x,y,:)) - meanTissue);
        end
    end
    
     figure(seriesNum*10 + 2*sliceNum)
     subplot(2,2,4), cla;
     imagesc(mydistimage), title('Manhattan distance to mean tissue curve'), colorbar;
    
    myslopefinder = -diff(fspecial('gaussian',[(st+1),1],7));
    pixelwiseUpslope = zeros(sx,sy);
    for x=1:sx
        for y=1:sy
            if(prelimMyo(x,y) > 0)
                pixelwiseUpslope(x,y) = max(filter2(myslopefinder,squeeze(cinemri1(x,y,:))));
            end
        end
    end
    figure(seriesNum*10 + 2*sliceNum+1), subplot(2,2,1)
    imagesc(pixelwiseUpslope), title('Pixelwise ratio of tissue upslope/AIF''s'); colorbar
    
    
    
    tissueDistmean = median(mydistimage(prelimMyo>0));
    tissueDiststd = std(mydistimage(prelimMyo>0));
    
    upslopemean = median(pixelwiseUpslope(prelimMyo>0));
    upslopestd = std(pixelwiseUpslope(prelimMyo>0));
    
    if(exist('pixelwisevariance'))
        variancemean = median(pixelwisevariance(prelimMyo>0));
        variancestd = std(pixelwisevariance(prelimMyo>0));
    end
    
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
    distMask = (mydistimage < (tissueDistmean + 1.5*tissueDiststd)) .* (mydistimage > (tissueDistmean - 1.5*tissueDiststd));
    b1Mask = (b1 < (b1mean + 1.5*b1std)) .* (b1 > (b1mean - 1.5*b1std));
    b2Mask = (b2 < (b2mean + 1.5*b2std)) .* (b2 > (b2mean - 1.5*b2std));
    b3Mask = (b3 < (b3mean + 1.5*b3std)) .* (b3 > (b3mean - 1.5*b3std));
    
    if(exist('pixelwisevariance'))
        varianceMask = (b3 < (b3mean + 1.5*b3std)) .* (b3 > (b3mean - 1.5*b3std));
    end
    upslopeMask = (pixelwiseUpslope < (upslopemean + 1.5*upslopestd)) .* (pixelwiseUpslope > (upslopemean - 1.5*upslopestd));
    
    delayWeigth = exp(-(mydelays -delaymean).^2 /(2*(1.5*delaystd)^2));
    distWeigth = exp(-(mydistimage -tissueDistmean).^2 /(2*(1.5*tissueDiststd)^2));
    b1Weigth = exp(-(b1 -b1mean).^2 /(2*(2*b1std)^2));
    b2Weigth = exp(-(b2 -b2mean).^2 /(2*(2*b2std)^2));
    b3Weigth = exp(-(b3 -b2mean).^2 /(2*(2*b3std)^2));
    
    if(exist('pixelwisevariance'))
        varianceWeight = exp(-(pixelwisevariance -variancemean).^2 /(2*(2*variancestd)^2));
    end
    upslopeWeight = exp(-(pixelwiseUpslope -upslopemean).^2 /(2*(2*upslopestd)^2));
    
    if(mean(delayMask(prelimMyo>0)) > .5)
        %prelimMyo = prelimMyo .* delayMask;
        prelimMyoWeight = prelimMyoWeight .* delayWeigth;
    end
    if(mean(distMask(prelimMyo>0)) > .5)
        %prelimMyo = prelimMyo .* distWeigth;
        prelimMyoWeight = prelimMyoWeight .* distWeigth;
    end
    if(mean(b1Mask(prelimMyo>0)) > .5)
        %prelimMyo = prelimMyo .* b1Mask;
        prelimMyoWeight = prelimMyoWeight .* b1Weigth;
    end
    if(mean(b2Mask(prelimMyo>0)) > .5)
        %prelimMyo = prelimMyo .* b2Mask;
        prelimMyoWeight = prelimMyoWeight .* b2Weigth;
    end
    if(mean(b3Mask(prelimMyo>0)) > .5)
        %prelimMyo = prelimMyo .* b3Mask;
        prelimMyoWeight = prelimMyoWeight .* b3Weigth;
    end
    if(mean(masterMask(prelimMyo>0)) > .5)
        %prelimMyo = prelimMyo .* masterMask;
        %prelimMyoWeight = prelimMyoWeight .* masterMask;
    end
    
    if(exist('pixelwisevariance'))
        if(mean(varianceMask(prelimMyo>0)) > .5)
            %prelimMyo = prelimMyo .* varianceMask;
            prelimMyoWeight = prelimMyoWeight .* varianceWeight;
        end
    end
    if(mean(upslopeMask(prelimMyo>0)) > .5)
        %prelimMyo = prelimMyo .* upslopeMask;
        prelimMyoWeight = prelimMyoWeight .* upslopeWeight;
    end
    mygaussian = fspecial('gaussian',[sx,sy],1);
    temp = filter2(mygaussian,prelimMyoWeight);
    temp2 = temp(prelimMyo>0); %#ok<NASGU>
    myo = temp>myothreshold;
    if(any(isnan(myo(:))) || any(isnan(prelimMyoWeight(:))))
        disp('Auto-segmentation falied.  Please run stage 3.11 for manual creation of contours');
        disp(['now Using default segmentation on series=' num2str(seriesNum) ' slice=' num2str(sliceNum)]);
                
                 theta = 0:.1:(2*pi);  % this part changed from just a return 9/19/11, EVRD
                AIF = horzcat((3*sin(theta)+sx/2)',(3*cos(theta)+sy/2)');
                epi = horzcat((10*sin(theta)+sx/2)',(10*cos(theta)+sy/2)');
                endo = horzcat((7*sin(theta)+sx/2)',(7*cos(theta)+sy/2)');
                angle = pi;
                save(['Output/blood_polyCoords.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'AIF','-ascii','-double');
                save(['Output/epi_polyCoords.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'epi','-ascii','-double');
                save(['Output/endo_polyCoords.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'endo','-ascii','-double');
                save(['Output/Roi_start_angle.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'angle','-ascii','-double');    
        return;
    end
    
    figure(seriesNum*10 + 2*sliceNum+1), subplot(2,2,2), imagesc(prelimMyoWeight,[0 median(prelimMyoWeight(prelimMyoWeight>0))+3*std(prelimMyoWeight(prelimMyoWeight>0))]), colorbar, title(['Weighted myocardium, threshold@ ' num2str(myothreshold)]);
    
    
    %smooth the myocardium
    myo = filter2(fspecial('gaussian',[sx,sy],1),myo)>.5;
    
    
    
    figure(seriesNum*10 + 2*sliceNum+1),subplot(2,2,3)
    myrgb = zeros(sx,sy,3);
    temp = cinemri1(:,:,22);
    temp = temp - min(temp(:));
    temp = temp/max(temp(:));
    myrgb(:,:,2) = temp;
    myrgb(:,:,3) = temp;
    temp(myo > 0) = min(1,temp(myo > 0)*2);
    myrgb(:,:,1) = temp;
    imagesc(myrgb);
    title('Myocardium')
    
    [r,c] = find(myo>0);
    [theta,radius] = cart2pol(r-meanLV(1),c- meanLV(2));
    if(isempty(theta))
        disp(['Myocardium was too hard to find.  Please perform manual segmentation on series=' num2str(seriesNum) ' slice=' num2str(sliceNum)]);
        disp(['Myocardium was too hard to find.  Using default segmentation on series=' num2str(seriesNum) ' slice=' num2str(sliceNum)]);
                
                 theta = 0:.1:(2*pi);  % this part changed from just a return 9/19/11, EVRD
                 AIF = horzcat((3*sin(theta)+sx/2)',(3*cos(theta)+sy/2)');
                epi = horzcat((10*sin(theta)+sx/2)',(10*cos(theta)+sy/2)');
                endo = horzcat((7*sin(theta)+sx/2)',(7*cos(theta)+sy/2)');
                angle = pi;
                save(['Output/blood_polyCoords.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'AIF','-ascii','-double');
                save(['Output/epi_polyCoords.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'epi','-ascii','-double');
                save(['Output/endo_polyCoords.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'endo','-ascii','-double');
                save(['Output/Roi_start_angle.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'angle','-ascii','-double');  
        return;
    end
  connectedMyo = EnsureConectivity(myo,meanLV);
    if(isempty(connectedMyo))
        disp(['Myocardium was too hard to find.  Please perform manual segmentation on series=' num2str(seriesNum) ' slice=' num2str(sliceNum)]);
        disp(['Myocardium was too hard to find.  Using default segmentation on series=' num2str(seriesNum) ' slice=' num2str(sliceNum)]);
                
                 theta = 0:.1:(2*pi);  % this part changed from just a return 9/19/11, EVRD
                 AIF = horzcat((3*sin(theta)+sx/2)',(3*cos(theta)+sy/2)');
                epi = horzcat((10*sin(theta)+sx/2)',(10*cos(theta)+sy/2)');
                endo = horzcat((7*sin(theta)+sx/2)',(7*cos(theta)+sy/2)');
                angle = pi;
                save(['Output/blood_polyCoords.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'AIF','-ascii','-double');
                save(['Output/epi_polyCoords.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'epi','-ascii','-double');
                save(['Output/endo_polyCoords.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'endo','-ascii','-double');
                save(['Output/Roi_start_angle.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'angle','-ascii','-double');    
        return;
    end
%     figure(seriesNum*10 + 2*sliceNum)
%     subplot(2,2,4), cla;
%     temp = parameters(:,:,1) - parameters(:,:,2).*parameters(:,:,3);
%     medk = median(temp(connectedMyo>0));
%     stdk = std(temp(connectedMyo>0));
%     imagesc(temp,[(medk-2*stdk) (medk+2*stdk)]), colormap jet; colorbar;title('Linear k-trans');
    
    
    %find optimal AIF contour
    LVCircle = [];
    maxVal = -1;
    AIFcurve = [];
    for x=1:sx
        for y=1:sy
            if(LVMask(x,y) > 0)
                theta = 0:.1:(2*pi);
                bloodX = 3*sin(theta) + x;
                bloodY = 3*cos(theta) + y;
                tempCircle = roipoly(LVMask,bloodX,bloodY);
                AIF = zeros(st,1);
                for t=1:st
                    temp = cinemri1(:,:,t);
                    AIF(t) = mean(temp(tempCircle>0));
                end
                if(max(AIF) > maxVal)
                    AIFcurve = AIF;
                    maxVal = max(AIF);
                    LVCircle = tempCircle;
                end
            end
        end
    end
    
    
    
    
    myslopefinder = -diff(fspecial('gaussian',[(st+1),1],7));
    AIF_upslope = max(filter2(myslopefinder,squeeze(AIFcurve)));
    pixelwiseUpslope = zeros(sx,sy);
    for x=1:sx
        for y=1:sy
            if(connectedMyo(x,y) > 0)
                pixelwiseUpslope(x,y) = max(filter2(myslopefinder,squeeze(cinemri1(x,y,:))));
            end
        end
    end
    pixelwiseUpslope = pixelwiseUpslope ./ AIF_upslope;
    figure(seriesNum*10 + 2*sliceNum+1), subplot(2,2,1)
    imagesc(pixelwiseUpslope), title('Pixelwise ratio of tissue upslope/AIF''s'); colorbar
    
    for i=1:1  %EVRD  bug SOMETIMES with 2, not sure need it? dilates more? 12/20/11
        %endocardium
        inner = imfill(~connectedMyo,round(meanLV));

        se = strel('disk',1);
        inner = imdilate(inner,se);
        [r,c] = find(inner>0);
        [~, minri] = min(r);
        boundary = bwtraceboundary(inner,[r(minri) c(minri)],'N');
        [theta,r] = cart2pol(boundary(:,1) -meanLV(1), boundary(:,2) -meanLV(2));
        r = r * 1.01;
        [x,y] = pol2cart(theta,r);
        x = x + meanLV(1); y = y + meanLV(2);
        endo = horzcat(x,y);
        endoLength=size(endo,1);
        if endoLength<30  %EVRD 12/20/11 to get around crash
            extendedEndo = vertcat(endo((end-1):end,:),endo,endo(1:endoLength,:));
        else
            extendedEndo = vertcat(endo((end-29):end,:),endo,endo(1:30,:));
        end
        clear smoothedExtendedEndo
        smoothedExtendedEndo(:,1) = conv(fspecial('gaussian',[length(extendedEndo) 1],3),extendedEndo(:,1),'same');
        smoothedExtendedEndo(:,2) = conv(fspecial('gaussian',[length(extendedEndo) 1],3),extendedEndo(:,2),'same');
        smoothedEndo = smoothedExtendedEndo(31:(end-30),:); %#ok<%#ok<MSNU> NASGU>
        endo = smoothedEndo;

        %epicardium
        endoMask = connectedMyo + inner;

        se = strel('disk',5);
        endoMask = imerode(endoMask,se);
        endoMask = imdilate(endoMask,se);
        endoMask = imdilate(endoMask,se);
        endoMask = imerode(endoMask,se);
        se = strel('disk',1);
        endoMask = imerode(endoMask,se);
        %endoMask = bwmorph(bwmorph(endoMask,'erode',5),'dilate',5);
        %endoMask = bwmorph(bwmorph(endoMask,'dilate',5),'erode',5);
        [r,c] = find(endoMask>0);
        [~, minri] = min(r);
        boundary = bwtraceboundary(endoMask,[r(minri) c(minri)],'N');
        [theta,r] = cart2pol(boundary(:,1) -meanLV(1), boundary(:,2) -meanLV(2));
        r = r * 1.01;
        [x,y] = pol2cart(theta,r);
        x = x + meanLV(1); y = y + meanLV(2);
        epi = horzcat(x,y);
        extendedepi = vertcat(epi((end-29):end,:),epi,epi(1:30,:));
        clear smoothedExtendedepi
        smoothedExtendedepi(:,1) = conv(fspecial('gaussian',[length(extendedepi) 1],1),extendedepi(:,1),'same');
        smoothedExtendedepi(:,2) = conv(fspecial('gaussian',[length(extendedepi) 1],1),extendedepi(:,2),'same');
        smoothedepi = smoothedExtendedepi(31:(end-30),:);
        epi = smoothedepi;
        if(i==1)
            figure(1);
            inner = roipoly(cinemri1(:,:,22),endo(:,2),endo(:,1));
            outter = roipoly(cinemri1(:,:,22),epi(:,2),epi(:,1));
            myo = outter - inner;
            connectedMyo = EnsureConectivity(myo,meanLV);
            if(isempty(connectedMyo))
                disp(['Myocardium was too hard to find.  Using default segmentation on series=' num2str(seriesNum) ' slice=' num2str(sliceNum)]);
                
                 theta = 0:.1:(2*pi);  % this part changed from just a return 9/19/11, EVRD
                AIF = horzcat((3*sin(theta)+sx/2)',(3*cos(theta)+sy/2)');
                epi = horzcat((10*sin(theta)+sx/2)',(10*cos(theta)+sy/2)');
                endo = horzcat((7*sin(theta)+sx/2)',(7*cos(theta)+sy/2)');
                angle = pi;
                save(['Output/blood_polyCoords.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'AIF','-ascii','-double');
                save(['Output/epi_polyCoords.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'epi','-ascii','-double');
                save(['Output/endo_polyCoords.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'endo','-ascii','-double');
                save(['Output/Roi_start_angle.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'angle','-ascii','-double');    
            %LVCircle = roipoly(LVMask,3*sin(theta)+sx/2)',(3*cos(theta)+sy/2)')
                return;
            end
        end
    end
    %blood pool
    [r,c] = find(LVCircle>0);
    [~, minri] = min(r);
    boundary = bwtraceboundary(LVCircle,[r(minri) c(minri)],'N');
    [theta,r] = cart2pol(boundary(:,1) -meanLV(1), boundary(:,2) -meanLV(2));
    r = r * 1.01;
    [x,y] = pol2cart(theta,r);
    x = x + meanLV(1); y = y + meanLV(2);
    blood = horzcat(x,y);
    
    if(~UseExistingAsPreliminary)
        %find the angle to the LV RV insertion point
        [angle,~] = cart2pol(topRV(1) - meanLV(1),topRV(2) - meanLV(2));
        %figure(seriesNum*10 + 2*sliceNum), subplot(2,2,1)
        %[x y] = pol2cart([(angle) (angle)],[0 norm([topRV(1) - meanLV(1),topRV(2) - meanLV(2)])]);
        %line(x+meanLV(1),y+meanLV(2));
    else
        angle = load(['Output/Roi_start_angle.study' int2str(seriesNum) '.slice' int2str(sliceNum)]);
    end
    if(isempty(epi) || isempty(endo) || isempty(angle))
        disp('Contours were not created, hmm- should have made default circles');
         
                 theta = 0:.1:(2*pi);  % this part changed from just a return 9/19/11, EVRD
                AIF = horzcat((3*sin(theta)+sx/2)',(3*cos(theta)+sy/2)');
                epi = horzcat((10*sin(theta)+sx/2)',(10*cos(theta)+sy/2)');
                endo = horzcat((7*sin(theta)+sx/2)',(7*cos(theta)+sy/2)');
                angle = pi;
                save(['Output/blood_polyCoords.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'AIF','-ascii','-double');
                save(['Output/epi_polyCoords.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'epi','-ascii','-double');
                save(['Output/endo_polyCoords.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'endo','-ascii','-double');
                save(['Output/Roi_start_angle.study' num2str(seriesNum) '.slice' num2str(sliceNum)],'angle','-ascii','-double');    
            
       % keyboard;
    end