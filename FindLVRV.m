
function [RV, LV] = findLVRV(cinemri1,showme)
if(~exist('showme'))
    showme = 0;
end


[sx sy st] = size(cinemri1);
smoothedcinemri1 = zeros([sx sy st]);
h = fspecial('gaussian',size([sx sy]), 1);
for t=1:st
    smoothedcinemri1(:,:,t) = filter2(h,cinemri1(:,:,t));
end
mymax = squeeze(max(smoothedcinemri1,[],3));
smoothed = filter2(fspecial('gaussian',size(mymax), 4),mymax);
BW = imregionalmax(smoothed);
if(showme)
    figure(1)
    clf;
    subplot(2,2,1)
    imagesc(mymax), colormap gray
    hold on
    subplot(2,2,2)
    imagesc(smoothed), colormap gray, hold on
    subplot(2,2,3); hold on;
    subplot(2,2,4)
    surf(smoothed)
end
[x, y] = find(BW);

if(showme)
    figure(1),clf, subplot(2,2,1), imagesc(mymax),colormap gray,hold on;
    for i=1:length(x)
        plot(y(i), x(i),'o');
    end
    hold off
end

temporalSmoothing = -diff(fspecial('gaussian',[1 (st+10+1)],3));
maxUpslope = zeros(length(x),3);
for i=1:length(x)
    values(i,:) = squeeze(mean(mean(cinemri1(max(1,x(i)-5):min(sx,x(i)+5),max(1,y(i)-5):min(sy,y(i)+5),:),1),2));
    values(i,:) = values(i,:) - mean(values(i,1:5));
    mymaxes(i) = max(values(i,:));
    signal = [repmat(mean(values(i,1:5)),[1,5]) values(i,:) repmat(mean(values(i,(end-5):end)),[1,5])];
    smoothedCurve = filter2(temporalSmoothing,signal);
    signal = signal(6:(end-5));
    smoothedCurve = smoothedCurve(6:(end-5));
    [maxUpslope(i,1) maxUpslope(i,2)] = max(smoothedCurve);
    maxUpslope(i,3) = norm([(x(i) - sx/2) (y(i) - sy/2)]); 
end



score = exp(-(maxUpslope(:,3)).^2/(2*(sx/12)^2)).*(1-1./(1+(maxUpslope(:,1)/(max(maxUpslope(:,1)))/2).^2));
if(length(score) >2)
    %IDX = kmeans(score,3,'EmptyAction','singleton');
    IDX = score>(max(score)/2);
    ratio = 3;
    while(sum(IDX) < 2)
        IDX = score>(max(score)/ratio);
        ratio = ratio + 1;
    end
    [worst worsti] = min(score);
    NotWanted = (IDX == IDX(worsti));

    %threshold = max(score) / 5;
    %disp(num2str(find(score > threshold)));
    temp = x;
    temp(NotWanted) = [];
    NumberOfTries = 1;
    x(NotWanted) = [];
    y(NotWanted) = [];

    values(NotWanted,:) = [];
    maxUpslope(NotWanted,:) = [];
    score(NotWanted) = [];
end
IDX = kmeans(maxUpslope(:,2),2,'EmptyAction','singleton');

RVspot = find(IDX == 1);
weights = score(RVspot);
weights = weights / sum(weights);
tempx = x(RVspot);
tempy = y(RVspot);
bestPoint = find(maxUpslope(RVspot,1) == max(maxUpslope(RVspot,1)));
%RV = [tempx(bestPoint) tempy(bestPoint)];
RV = round([sum(tempx.*weights) sum(tempy.*weights)]);
LVspot = find(IDX == 2);
weights = score(LVspot);
weights = weights / sum(weights);
tempx = x(LVspot);
tempy = y(LVspot);
bestPoint = find(maxUpslope(LVspot,1) == max(maxUpslope(LVspot,1)));
%LV = [tempx(bestPoint) tempy(bestPoint)];
LV = round([sum(tempx.*weights) sum(tempy.*weights)]);
%switch LV/RV if LV has an earlier upslope
if(LV(2) < RV(2))
    temp = RV;
    RV = LV;
    LV = temp;
end

trueMax = max(maxUpslope(:,1));
colors = lines(length(x));
if(showme)
    for i=1:length(x)
        mycolor = colors(IDX(i),:);
%         myhsv = rgb2hsv(mycolor);
%         myhsv(2) = score(i)/bestScore;
%         mycolor = hsv2rgb(myhsv);
        subplot(2,2,1)
        plot(y(i), x(i),'o','Color',mycolor);
        plot(values(i,:), 'Color',mycolor);
    end
    
    subplot(2,2,1)
    hold off;
    subplot(2,2,3)
    hold off;
    
    figure(3)
    clf;
    
    subplot(2,2,1), 
    imagesc(cinemri1(:,:,22)), colormap gray
    hold on
    plot(RV(2), RV(1),'o');
    plot(LV(2), LV(1),'or');
    hold off
    
    subplot(2,2,2)
    hold on
    plot(squeeze(mean(mean(cinemri1((RV(1)-5):(RV(1)+5),(RV(2)-5):(RV(2)+5),:),1),2)));
    plot(squeeze(mean(mean(cinemri1((LV(1)-5):(LV(1)+5),(LV(2)-5):(LV(2)+5),:),1),2)),'r');
    hold off
    legend('RV','LV')
    
end
