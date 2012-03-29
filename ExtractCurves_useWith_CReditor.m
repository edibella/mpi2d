function curves = ExtractCurves_useWith_CReditor(cinemri1,blood,endo,epi,angle,numAzimuthalRegions,referenceFrame)
    %do a test on the angle to see if it's legacy(mod(angle,360)) or updated(mod(angle,2*pi))
    legacy = 0;
    if(angle > 2*pi)
        legacy = 1;
        angle = mod((angle-90)*pi/180,2*pi);
    end
    
    [sx,sy,st] = size(cinemri1);
    %use the contours to recreate a mask as that's what the user would be
    %doing to recreate our results
    outtermask = roipoly(cinemri1(:,:,22),epi(:,1),epi(:,2));
    innermask = roipoly(cinemri1(:,:,22),endo(:,1),endo(:,2));
    newMyo = (outtermask - innermask)>0;
    newLV = roipoly(cinemri1(:,:,22),blood(:,1),blood(:,2));
    %[r,c] = find(newLV>0);
    %newMeanLV = [mean(r) mean(c)];
    
    %use the myocardial mask to find the curves then the deltaSI curves
    [r,c] = find(outtermask>0);
    middle = [mean(r) mean(c)];
    [r,c] = find(newMyo>0);
    [theta,~] = cart2pol(r-middle(1), c - middle(2));
    theta = mod(theta - angle,2*pi);
    
    mycolors = reshape(lines(numAzimuthalRegions),1,numAzimuthalRegions,3);
    
    %subplot(2,2,3)   %EVRD
    %figure;
    myrgb = zeros(sx,sy,3);
    temp = mean(cinemri1(:,:,(referenceFrame-3):(referenceFrame+3)),3);
    temp = temp - min(temp(:));
    temp = temp/max(temp(:));
    myrgb(:,:,2) = temp;
    myrgb(:,:,3) = temp;
    myrgb(:,:,1) = temp;
    myregionalMask = zeros(sx,sy);
    labeled = floor(theta/(2*pi)*numAzimuthalRegions);
    %make it so that  region 3 and 4 are in the septal wall
    %(the region counter clockwise from the insertion angle is region 3, and the next clockwise region is region 4)
    if(~legacy)
        labeled = mod(labeled+2,6);
    end
    for i=1:length(r)
        myregionalMask(r(i),c(i)) = labeled(i)+1;
        myrgb(r(i),c(i),:) = myrgb(r(i),c(i),:).*mycolors(1,labeled(i)+1,:);
    end
    %imagesc(myrgb);
    
    for region=1:numAzimuthalRegions
        mymask = myregionalMask==region;
        [pointx,pointy] = find(mymask);
        pointx = mean(pointx);
        pointy = mean(pointy);
        text(pointy,pointx,num2str(region),'Color',[1 1 1]);
    end
    
    
    curves = zeros(st,numAzimuthalRegions);
    for region=1:numAzimuthalRegions
        regionalMask = myregionalMask == region;
        for t=1:st
            temp2 = cinemri1(:,:,t);
            curves(t,region) = mean(temp2(regionalMask>0));
        end
    end
    %other = load('/v/raid1/npack/Processing/P010710_nufft/Output/curves.study13.slice1.mat');
    %figure(2), subplot(2,2,1),plot(other.curves(2:end,:)'),subplot(2,2,2),plot(curves), subplot(2,2,3), plot(abs(other.curves(2:end,:)' - curves));
    %return;
    fifthHighestInAIF = zeros(st,1);
    for t=1:st
        temp = cinemri1(:,:,t);
        values = temp(newLV>0);
        values = sort(values);
        fifthHighestInAIF(t) = values(end-5);  % how similar to mean in the case P071709??
        meanInAIF(t) = mean(temp(newLV>0));
    end
    curves = vertcat(meanInAIF,curves');

    if(any(isnan(curves)))
        disp('Houston we have a curve problem');
        try
            if(exist('mpi2d.log'))
                fid = fopen('mpi2d.log','a'); 
            else
                fid = fopen('mpi2d.log','w'); 
            end
            fprintf(fid,'%s',['There were nans in the curves file typically implying that the contours led to invalid masks.  They were probably createdin stage 3.12./n']); 
            fclose(fid);
        catch ME %#ok<NASGU>

        end
        %keyboard;
    end
    %subplot(2,2,4), cla,hold on, plot(curves(1,:)','k');  %EVRD 12/15/10
%     figure; plot(curves(1,:)','k');
%     hold on
%     for region=1:numAzimuthalRegions
%         plot(curves(region+1,:)','Color',squeeze(mycolors(1,region,:)));
%     end
%     title('Tissue+AIF curves');
%     legend('AIF','1','2','3','4','5','6')
    
    return;
    
    