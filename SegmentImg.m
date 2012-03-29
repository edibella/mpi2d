function segmented = SegmentImg(ktransImg, DeltaTImg,cine)

[s1 s2] = size(ktransImg);
sigmaIntensity = .2;
intensities = ktransImg(:);
[Xi Yi] = meshgrid(intensities,intensities);
wintensity1 = exp(-(Xi - Yi).^2./(2*sigmaIntensity^2));


[s1 s2] = size(DeltaTImg);
sigmaIntensity = 2;
intensities = DeltaTImg(:);
[Xi Yi] = meshgrid(intensities,intensities);
wintensity2 = exp(-(Xi - Yi).^2./(2*sigmaIntensity^2));


affinityMatrix = wintensity1.*wintensity2;

%normalization
% D = sum(affinityMatrix,2);
% D = diag(D);
% D2 = D^(-1/2);
% affinityMatrix = D2*(affinityMatrix)*D2;
% if(isnan(sum(affinityMatrix(:))) || isinf(sum(affinityMatrix(:))))
%     disp('We''ve got a nan');
%     return;
% end

subplot(2,3,2), imagesc(affinityMatrix), colormap gray

disp('Computing eigenvector...');
pause(.05)


[V,D] = eig(affinityMatrix);
D = diag(D);
[Dsorted IX] = sort(D,'descend');
subplot(2,3,3), plot(Dsorted), title('Eigenvalues');
runningimg = zeros(s1, s2,3);
segmented = zeros(s1,s2);
subplot(2,3,6);
plot(0);
hold on;
runningSegmentaiton = [];
clear randcolor;
segment = 1;
randcolor = [];
while (sum(abs(V(:))) > 0)
    [dmax index]= max(D);
    D(index) = 0;
    vector = V(:,index);
    temp = abs(vector);
    goodPoints = temp > mean(temp(temp>0))/1.1;
    if(sum(goodPoints) < 20)
        goodPoints = ~(sum(runningSegmentaiton,2)>0);
    end
    runningSegmentaiton(:,end+1) = goodPoints;
    V(goodPoints,:) = 0;
    if(sum(goodPoints) == 0) 
        break;
    end
    randcolor(end+1,:) = 1-(rand(1,3).^2);
    x = find(goodPoints > 0);
    y = temp(x);
    subplot(2,3,6);
    plot(x,y,'x','Color',randcolor(end,:),'MarkerSize',3);
    title('Colored Eigenvectors, matches segmentation');
    
    mask = reshape(goodPoints,[s1 s2]);
    [x y] = find(mask);
    for i=1:length(x)
        runningimg(x(i),y(i),1) = randcolor(end,1);
        runningimg(x(i),y(i),2) = randcolor(end,2);
        runningimg(x(i),y(i),3) = randcolor(end,3);
        segmented(x(i),y(i)) = segment;
    end
    subplot(2,3,4),imagesc(runningimg), colormap gray;
    title('Segmentation');
    segment = segment + 1;

end

%reorder afinity matrix
%reorder horizontally
storage = [];
for i=1:size(runningSegmentaiton,2)
    index = size(storage,1);
    storage((index+1):(index + sum(runningSegmentaiton(:,i))),:) = affinityMatrix(runningSegmentaiton(:,i)>0,:);
end
affinityMatrix = storage;

%reorder vertically
storage = [];
for i=1:size(runningSegmentaiton,2)
    index = size(storage,2);
    storage(:,(index+1):(index + sum(runningSegmentaiton(:,i)))) = affinityMatrix(:,runningSegmentaiton(:,i)>0);
end
affinityMatrix = storage;
subplot(2,3,5),imagesc(affinityMatrix), colormap gray;
title('Sorted Afinity Matrix');

disp('Done');