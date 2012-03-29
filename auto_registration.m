function offset = auto_registration(matfilename_or_cinemri,ranget,UserInterventionBoolean) 
mywaitbar = waitbar(0,'Finding Edges');
if(size(matfilename_or_cinemri,3) == 1)
    disp(['Fixing: ' matfilename_or_cinemri]);
    load(matfilename_or_cinemri);
    if(exist('imgrr'))
        originalImage = imgrr(:,:,ranget);
    elseif(exist('cinemri1'))
        originalImage = cinemri1(:,:,ranget);
    elseif(exist('cinemri'))
        originalImage = cinemri(:,:,ranget);
    end
else
    [s1 s2 s3] = size(matfilename_or_cinemri);
    disp(['Fixing: a ' num2str(s1) 'x' num2str(s2) 'x' num2str(s3) ' movie']);
    originalImage = matfilename_or_cinemri;
end

maxa = max(originalImage(:));
[s1 s2 tmax] = size(originalImage);

% figure(2);
% imagesc(a(:,:,22));
% colormap gray;
% disp('Please outline a box around something you want to remain still');
% [BW, xp, yp] = roipoly




edges = zeros(s1,s2,tmax);
tedges = zeros(s1,s2);
originalImage = originalImage-mean(originalImage(:));

for t=1:tmax
    cropped2(:,:,t) = wiener2(originalImage(:,:,t));
end
mincrop = min(cropped2(:));
maxcrop = max(cropped2(:));
tedges = zeros(size(cropped2,1),size(cropped2,2));
for t=1:tmax
    mywaitbar = waitbar(t/tmax);
    % figure(6);
%     subplot(1,2,1);
%     imagesc(cropped2(:,:,t),[mincrop maxcrop]);
%     colormap gray;
%     subplot(1,2,2);
%     imagesc(edges(:,:,t));
%     colormap gray;
%     pause(.1);
    edges(:,:,t) =  edge(cropped2(:,:,t),'canny',[.03 .3]);
    tedges = tedges + edges(:,:,t);
    [fx fy] = gradient(cropped2(:,:,t)); 
    gradients(:,:,t) = sqrt(fx.*fx + fy.*fy);

    %cropped(:,:,t) = gradients(:,:,t);
end
% h = figure(2);
% close(h);
% h = figure(2);
% imagesc(tedges>5);
% colormap gray;
% disp('Please outline a box around something you want to remain still');
% [BW, xp, yp] = roipoly;
% close(h);

BW = roipoly(tedges>5,size(tedges,1)/2*cos(0:.05:2*pi)+size(tedges,1)/2,size(tedges,2)/2*sin(0:.05:2*pi)+size(tedges,2)/2);


%figure(5),imagesc(tedges),colormap gray;pause;


h = fspecial('gaussian',[7 7],5);
h2 = fspecial('gaussian',[15 15],5);
clear temp;
%do the first and last frames
for t = [1 2 tmax tmax-1]
    %stuffToCorrelate(:,:,t) = filter2(h2,BW).*cropped2(:,:,t) .*(filter2(h,edges(:,:,t))>0);
    stuffToCorrelate(:,:,t) = BW.*edges(:,:,t);
    %temp = stuffToCorrelate(:,:,t);
    %temp = mean(temp(temp>0));
    %stuffToCorrelate(:,:,t) = stuffToCorrelate(:,:,t) + (~BW).*temp;
    %figure(10),imagesc(stuffToCorrelate(:,:,t)),colormap gray, pause(.1);
    
end
%blur each frame's edge image
for t=3:tmax-2
    %stuffToCorrelate(:,:,t) = filter2(h2,BW).*cropped2(:,:,t).*(filter2(h,edges(:,:,t))>0);
    stuffToCorrelate(:,:,t) = BW.*edges(:,:,t);
    %stuffToCorrelate(:,:,t) = filter2(h2,  stuffToCorrelate(:,:,t));
    %figure(10),imagesc(stuffToCorrelate(:,:,t)),colormap gray, pause(.1);
end

%composite this with adjacent frames into temporary matrix
% stuffToCorrelate2 = zeros(size(stuffToCorrelate));
% for t=3:tmax-2
%     stuffToCorrelate2(:,:,t) = stuffToCorrelate(:,:,t).*cropped2(:,:,t) + .3*stuffToCorrelate(:,:,t-1)+ .3*stuffToCorrelate(:,:,t+1) + .1*stuffToCorrelate(:,:,t-2)+ .1*stuffToCorrelate(:,:,t+2);
%     stuffToCorrelate2(:,:,t) = stuffToCorrelate2(:,:,t)/(max(max(stuffToCorrelate2(:,:,t))));
%     temp(:,:,t) = stuffToCorrelate2(:,:,t).*edges(:,:,t);
%     stuffToCorrelate2(:,:,t) = temp(:,:,t);% + (1-stuffToCorrelate2(:,:,t))*mean(temp(:));
%     figure(10),imagesc(stuffToCorrelate2(:,:,t)),colormap gray,pause(.1);
% end
% clear temp;
% %store it back to original matrix
% for t=3:tmax-2
%     stuffToCorrelate(:,:,t) = stuffToCorrelate2(:,:,t);
% end



[s1 s2 s3] = size(stuffToCorrelate);
corr_offset = zeros(tmax,2);

c = normxcorr2(stuffToCorrelate(:,:,1),stuffToCorrelate(:,:,1+1));
h = fspecial('gaussian',[size(c,1),size(c,2)],size(c,1)/6);
howFarBack = min(7,s3-2);

%create triangle from back to front and align up the last howFarBack frames
counter = 1;
close(mywaitbar);
mywaitbar = waitbar(0,'Registering Image: Backward Alignment');
clear max_cs
shift = zeros(tmax,2);
for t=tmax-1:-1:tmax-howFarBack-1
    waitbar((tmax-t)/tmax);
    max_cs = zeros(counter,1);
    for i=1:counter
        c = normxcorr2(stuffToCorrelate(:,:,t),stuffToCorrelate(:,:,t+i));
        temp = c.*h;
        [max_c, imax] = max(abs(temp(:)));
        max_cs(i) = max_c;
        [ypeak, xpeak] = ind2sub(size(c),imax(1));
        corr_offset(t,:,i) = [(ypeak-s1) (xpeak-s2)];
    end
    if(counter == 1)
        shift(t,:) = corr_offset(t,:,:);
        %standardDeviations(t) = 0;
    else
        shift(t,:) = floor(median(squeeze(corr_offset(t,:,:))'));
        %standardDeviations(t) = std(norm(squeeze(corr_offset(t,:,:))'));
        %standardDeviations(t) = std(max_cs);
    end
    cropped2(:,:,t) = circshift(cropped2(:,:,t),shift(t,:));
    gradients(:,:,t) = circshift(gradients(:,:,t),shift(t,:));
    stuffToCorrelate(:,:,t) = circshift(stuffToCorrelate(:,:,t),shift(t,:));
    counter = counter+1;
end



%now that the last howFarBack are aligned we can use them to align in a
%window going to the first frame

clear max_cs
max_cs = zeros(howFarBack,1);
corr_offset2 = zeros(tmax,2,howFarBack);
standardDeviation = std((shift((size(shift,1)-howFarBack):end,1)).^2 + (shift((size(shift,1)-howFarBack):end,2)).^2);
if(standardDeviation == 0) standardDeviation = size(c,2)/10; end;
hh = fspecial('gaussian',[size(c,1),size(c,2)],8*standardDeviation);

for t=tmax-howFarBack:-1:1
    mywaitbar = waitbar((tmax-t)/tmax);
    for i=1:howFarBack
        try
            c = normxcorr2(stuffToCorrelate(:,:,t),stuffToCorrelate(:,:,t+i));
            blending_factor = 1/(1+((t/tmax)/.3)^2);%at about 20% of the way to the first frame transition to focusing more on only considering shifts like what we've already seen
            temp = c.*(h*(1-blending_factor) + blending_factor*hh);

    %         if(t < 30)
    %             figure(3);
    %             title([num2str(t) 'vs. ' num2str(i)]);
    %             subplot(2,2,1);
    %             surf(temp);
    %             subplot(2,2,3);
    %             imagesc(stuffToCorrelate(:,:,t));
    %             colormap gray
    %             title([num2str(t)]);
    %             subplot(2,2,2);
    %             surf((h*(1-blending_factor) + blending_factor*hh));
    %             subplot(2,2,4);
    %             imagesc(stuffToCorrelate(:,:,t+i));
    %             colormap gray
    %             title([num2str(i)]);
    %             pause(.05);
    %         end

            [max_c, imax] = max(abs(temp(:)));
            max_cs(i) = max_c;
            [ypeak, xpeak] = ind2sub(size(c),imax(1));
            corr_offset2(t,:,i) = [(ypeak-s1) (xpeak-s2)];
        catch
            %The stuffToCorrelate matrix had nothing in it/had nothing as
            %reference.  I'll simply assume there's zero shift
            corr_offset2(t,:,i) = [0 0];
        end
        
    end
    
    %standardDeviations(t) = std(norm(squeeze(corr_offset2(t,:,:))'));
    %standardDeviations(t) = std(max_cs);
    shift(t,:) = floor(median(squeeze(corr_offset2(t,:,:))'));
    cropped2(:,:,t) = circshift(cropped2(:,:,t),shift(t,:));
    gradients(:,:,t) = circshift(gradients(:,:,t),shift(t,:));
    stuffToCorrelate(:,:,t) = circshift(stuffToCorrelate(:,:,t),shift(t,:));
    
    %recalculate the standard deviation and hh weighting function
    standardDeviation = std((shift((size(shift,1)-howFarBack):end,1)).^2 + (shift((size(shift,1)-howFarBack):end,2)).^2);
    if(standardDeviation == 0) standardDeviation = size(c,2)/10; end;
    hh = fspecial('gaussian',[size(c,1),size(c,2)],8*standardDeviation);
end
close(mywaitbar);

% %repeat but going forward this time
% standardDeviations = zeros(tmax,1);
% mywaitbar = waitbar(0,'Registering Image: Forward Pass');
% counter = 1;
% clear corr_offset
% disp('Creating header for forward pass');
% for t=2:howFarBack
%     mywaitbar = waitbar(t/tmax);
%     for i=1:counter
%         c = normxcorr2(stuffToCorrelate(:,:,t),stuffToCorrelate(:,:,t-i));
%         temp = c.*h;
%         [max_c, imax] = max(abs(temp(:)));
%         [ypeak, xpeak] = ind2sub(size(c),imax(1));
%         corr_offset(t,:,i) = [(ypeak-s1) (xpeak-s2)];
%     end
%     if(counter == 1)
%         shift = corr_offset(t,:,:);
%         %standardDeviations(t) = 0;
%     else
%         shift = floor(median(squeeze(corr_offset(t,:,:))'));
%         %standardDeviations(t) = std(norm(squeeze(corr_offset(t,:,:))'));
%         %standardDeviations(t) = std(max_cs);
%     end
%     
%     cropped2(:,:,t) = circshift(cropped2(:,:,t),shift);
%     gradients(:,:,t) = circshift(gradients(:,:,t),shift);
%     stuffToCorrelate(:,:,t) = circshift(stuffToCorrelate(:,:,t),shift);
%     counter = counter+1;
% end
% 
% disp('Doing forward pass');
% max_cs = zeros(howFarBack,1);
% corr_offset2 = zeros(tmax,2,howFarBack);
% for t=howFarBack+1:tmax
%     mywaitbar = waitbar(t/tmax);
%     for i=1:howFarBack
%         c = normxcorr2(stuffToCorrelate(:,:,t),stuffToCorrelate(:,:,t-i));
%         temp = c.*h;
%         [max_c, imax] = max(abs(temp(:)));
%         max_cs(i) = max_c;
%         [ypeak, xpeak] = ind2sub(size(c),imax(1));
%         corr_offset2(t,:,i) = [(ypeak-s1) (xpeak-s2)];
%     end
%    
%     
%     %standardDeviations(t) = std(max_cs);
%     shift = floor(median(squeeze(corr_offset2(t,:,:))'));
%     cropped2(:,:,t) = circshift(cropped2(:,:,t),shift);
%     gradients(:,:,t) = circshift(gradients(:,:,t),shift);
%     stuffToCorrelate(:,:,t) = circshift(stuffToCorrelate(:,:,t),shift);
% end
%close(mywaitbar);



tedges = zeros(size(cropped2,1),size(cropped2,2));
for t=1:tmax
    edges(:,:,t) =  edge(cropped2(:,:,t),'canny',[.07 .3]);
    tedges = tedges + edges(:,:,t);
end


if(UserInterventionBoolean)
    %reconstruct the shifts necessary to get where we are then subtract the DC
    %component
    clear global values
    global values
    values.tedges = tedges;
    values.BW = BW;
    values.threshold = 15;
    values.mask = (values.tedges>values.threshold ).*values.BW;
    disp('PageDown/PageUp to go forward/backward in time.');
    disp('Space and Enter also progress time');
    disp('Use arrow keys or wasd to move frame.');
    disp('Press escape when done');
    for t=1:tmax
        temp2 = cropped2(:,:,t);
        auto_contrast(:,:,t) = histeq(cropped2(:,:,t)/max(temp2(:)),200);
        c = normxcorr2(cropped2(:,:,t),originalImage(:,:,t));
        [max_c, imax] = max(abs(c(:)));
        max_cs(i) = max_c;
        [ypeak, xpeak] = ind2sub(size(c),imax(1));
        values.corr_offset(t,:) = [(ypeak-s1) (xpeak-s2)];
    end
    values.corr_offset= sqrt(values.corr_offset(:,1).^2 + values.corr_offset(:,2).^2);
    standardDeviations = values.corr_offset;
    %create a confidence interval
    values.alpha = .1;
    talphaover2 = tinv(1-values.alpha, size(values.corr_offset,3)-1);
    values.lowerbound = mean(values.corr_offset)-talphaover2*std(values.corr_offset)/sqrt(size(values.corr_offset,1));
    values.upperbound = mean(values.corr_offset)+talphaover2*std(values.corr_offset)/sqrt(size(values.corr_offset,1));
    outliers = values.corr_offset(values.corr_offset<values.lowerbound .* values.corr_offset<values.upperbound);
    if(size(outliers) > 0)
        %grab the first one and have the user look at it
        values.t = find(values.corr_offset==outliers(1));
    else
        values.t = 1;
    end

    values.done = 0;
    values.shift = [0 0];
    values.tmax = tmax;
    h = figure(1);
    close(h);
    h = figure(1);
    set(h, 'OuterPosition',[44   182  1093   807]);
    set(h,'KeyPressFcn',@printfig);
    range = [min(cropped2(:)) .7*max(cropped2(:))];
    while values.done == 0
        %t = values.t; 

        if(values.t > values.tmax) values.t = values.tmax; end
        if(values.t < 1) values.t = 1; end

        subplot(2,2,1);
        imagesc(max(values.mask*range(2),cropped2(:,:,values.t)));%,range);
        title('Original auto-Registered');
        colormap gray;axis image;
        subplot(2,2,2);
        imagesc(max(values.mask,auto_contrast(:,:,values.t)));
        title('Auto-contrast');
        colormap gray;axis image;
        subplot(2,2,4);
        imagesc(gradients(:,:,values.t));
        title('Gradient');
        colormap gray;axis image;
        subplot(2,2,3);
        plot(standardDeviations);
        title('Magnitude of the shifts');
        hold on;
        plot([values.t values.t],[0 max(standardDeviations)],'r');
        hold off;axis square;
        %imagesc(stuffToCorrelate(:,:,t),[min(stuffToCorrelate(:)) .6*max(stuffToCorrelate(:))]);
        %colormap gray;colorbar;
        if(norm(values.shift) > 0)
            cropped2(:,:,values.t) = circshift(cropped2(:,:,values.t),values.shift);
            auto_contrast(:,:,values.t) = circshift(auto_contrast(:,:,values.t),values.shift);
            gradients(:,:,values.t) = circshift(gradients(:,:,values.t),values.shift);
            stuffToCorrelate(:,:,values.t) = circshift(stuffToCorrelate(:,:,values.t),values.shift);
            standardDeviations(values.t) = 0;%user fixed it
            values.shift = [0 0];
        else
            values.t = values.t + 0;
        end
        if(values.t > tmax)
            values.t = 1;
        end
        pause(.05);
    end
    close(h);
end
%reconstruct the shifts necessary to get where we are
clear offset
for t=1:tmax
    c = normxcorr2(cropped2(:,:,t),originalImage(:,:,t));
    [max_c, imax] = max(abs(c(:)));
    max_cs(i) = max_c;
    [ypeak, xpeak] = ind2sub(size(c),imax(1));
    offset(t,:) = [(ypeak-s1) (xpeak-s2)];

end

%subtract the DC component
for t=1:tmax
    offset(t,:) = offset(t,:)- round(mean(offset));
end
%offset
end


function printfig(src,evnt)
    global values
    keyinput = evnt.Key;
    if strcmp(keyinput,'h') || strcmp(keyinput,'w') ||  strcmp(keyinput,'uparrow')
        values.shift = [-1 0];
    end
    if strcmp(keyinput,'j') || strcmp(keyinput,'d') ||  strcmp(keyinput,'rightarrow')
        values.shift = [0 1];
    end
    if strcmp(keyinput,'k') || strcmp(keyinput,'a') ||  strcmp(keyinput,'leftarrow')
        values.shift = [0 -1];
    end
    if strcmp(keyinput,'l') || strcmp(keyinput,'s') ||  strcmp(keyinput,'downarrow')
        values.shift = [1 0];
    end
    if strcmp(keyinput,'subtract')
        values.threshold = values.threshold+1;
        values.mask = (values.tedges>values.threshold ).*values.BW;
    end
    if strcmp(keyinput,'add')
        values.threshold = values.threshold-1;
        values.mask = (values.tedges>values.threshold ).*values.BW;
    end
    if strcmp(keyinput,'space') ||  strcmp(keyinput,'return') ||  strcmp(keyinput,'pagedown')
        values.t = values.t+1;
        if(values.t == values.tmax)
            values.t = 1;
        end
    end
    if strcmp(keyinput,'n')
        %are there any outliers
        outliers = values.corr_offset(values.corr_offset<values.lowerbound .* values.corr_offset<values.upperbound);
        if(size(outliers) > 0)
            %grab the first one and have the user look at it
            values.t = find(values.corr_offset==outliers(1));
        else
            values.alpha = values.alpha / 2;
            talphaover2 = tinv(1-values.alpha, size(values.corr_offset,3)-1);
            values.lowerbound = mean(values.corr_offset)-talphaover2*std(values.corr_offset)/sqrt(size(values.corr_offset,1));
            values.upperbound = mean(values.corr_offset)+talphaover2*std(values.corr_offset)/sqrt(size(values.corr_offset,1));

            disp(['Recreating shift bounds: ' num2str(values.lowerbound) ' - ' num2str(values.upperbound)]);
            if(size(outliers) > 0)
                %grab the first one and have the user look at it
                values.t = find(values.corr_offset==outliers(1));
            else
                values.t = values.t+1;
                if(values.t == values.tmax)
                    values.t = 1;
                end
            end
        end
    end
    if strcmp(keyinput,'backspace') ||  strcmp(keyinput,'pageup')
        values.t = values.t-1;
        if(values.t == 0)
            values.t = values.tmax;
        end
    end
    if strcmp(keyinput,'escape')
        values.done = 1;
    end
end

