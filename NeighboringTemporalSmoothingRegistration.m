function [shifts,updatedcine] = NeighboringTemporalSmoothingRegistration(cinemri1, rangex, rangey, Nitter, timestamps,show)
[sx sy st] = size(cinemri1);
saved = cinemri1;
[srangex srangey] = size(cinemri1(rangex, rangey, 1));
shifts = zeros(st,2);
cc = normxcorr2(cinemri1(rangex,rangey,2),cinemri1(rangex,rangey,1)); 
if(~exist('show'))
    show = 0;
end


b = (-log(1.5*st)/log(2) + Nitter)/(Nitter-1);
a = 1-b;
b = 2^b;
windowSizes = b*2.^(a*(1:Nitter));
shiftHistory = [];
for itter=1:Nitter
%     if(itter > Nitter/2)
%         oldshifts = shifts;
%         shifts = fixshifts(shifts, timestamps, 1);
%         
%         if(any(any(abs(shifts)>100)))
%             disp('Problem houston.  Shifts are gihugoumongeous');
%             keyboard;
%         end
%         for t=1:st
%             cinemri1(:,:,t) = circshift(saved(:,:,t),shifts(t,:));
%         end
%     end
    
    
    mask = fspecial('disk',20);%max(ceil(20*(Nitter-itter+1)/(Nitter+1)),10));
    temp = zeros(size(cc));
    range1 = (1+round(size(cc,1)/2-size(mask,1)/2));
    range1 = range1:(range1+size(mask,1)-1);
    range2 = (1+round(size(cc,2)/2-size(mask,2)/2));
    range2 = range2:(range2+size(mask,2)-1);
    temp(range1, range2) = mask;
    mask = ~temp;
    
    windowSize = round(windowSizes(itter));
    if(mod(windowSize,2) == 0)
        windowSize = windowSize + 1;
    end
    side = floor(windowSize/2);
    
    
    h = fspecial('gaussian',[1 windowSize],windowSizes(itter)/1.5);
    h(side+1) = 0;
    h(1:side) = 2*h(1:side);
    summ = sum(h);
    h = h/summ; %make it sum to 1
    h = reshape(h,1,1,windowSize);
    h = repmat(h,[srangex srangey]);
    
    padded = cat(3,repmat(mean(cinemri1(:,:,1:side),3),[1 1 side]), cinemri1, repmat(mean(cinemri1(:,:,(end-side):end),3),[1 1 side]));
    %order = randperm(st);
    order = 1:st;
    for t=order
        refrenceFrame = sum(padded(rangex, rangey, [((t-side):(t+side))+side]).*h,3);
        cc = normxcorr2(refrenceFrame,cinemri1(rangex,rangey,t)); 
        c = cc;  %c(mask) = 0;
        hh = zeros(size(c));
        width = max(1,round(min(t/2,20)));
        range1 = round(size(cc,1)/2-width):round(size(cc,1)/2+width);
        range2 = round(size(cc,2)/2-width):round(size(cc,2)/2+width);
        hh(range1,range2) = 1;
        %if(t < 10)
        %    figure(5),subplot(2,1,1),surf(cc),title([num2str(itter) ' - ' num2str(t)]), subplot(2,1,2), surf(cc.*hh)
        %end
        cc = cc .* hh;
        [max_cc, imax] = max(abs(cc(:)));
        [xpeak, ypeak] = ind2sub(size(cc),imax(1));
%         if(any(any(cc((xpeak-1):(xpeak+1), (ypeak-1):(ypeak+1)) > max_cc)))
%             disp([num2str(t) ' was not correlatable']);
%             continue;
%         end
        shiftAmount = -[ (xpeak-srangex) (ypeak-srangey) ];
        shiftHistory = [shiftHistory norm(shiftAmount)];
        if(norm(shiftAmount) > 15) 
            continue; 
        end
        shifts(t,:) = shifts(t,:) + shiftAmount;
        
        %if(t < 10)
        %    subplot(2,1,2), title([num2str(shifts(t,:))]); pause(.1)
        %end
        cinemri1(:,:,t) = circshift(cinemri1(:,:,t),shiftAmount);
    end
    
    %figure(1), subplot(2,2,4)
    %plot(shiftHistory,'.');
end
updatedcine = cinemri1;
if(show)
    for t=1:st
        figure(1)
        subplot(2,2,1), imagesc(saved(:,:,t)), colormap gray;
        subplot(2,2,2), imagesc(cinemri1(:,:,t)), colormap gray;
        subplot(2,2,3), imagesc(cinemri1(rangex,rangey,t)), colormap gray;
        pause(.1);
    end
end