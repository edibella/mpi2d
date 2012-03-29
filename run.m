%% make a filter
% fid = fopen('/home/mirl/bmatthew/bad files.txt');
%     
% tline = fgetl(fid);
% while ischar(tline)
%     A = sscanf(tline,'%f %f');
%     magnitude = sqrt(A(1)*A(1) + A(2)*A(2));
%     if(magnitude > 10)
%         %disp([num2str(magnitude) ' : ' directory '/' char(files(i).name)]);
%         flaggedFiles.badfiles(end+1) = {[ directory '/' char(files(i).name)]};
%         break;
%     end
%     tline = fgetl(fid);
% end
% fclose(fid);
load '/v/raid1/npack/MRIdata/Cardiac/Trio/P043009/ImageData/nufftP043009MID34_slice1_100iter_4coils_F1000_T900_S0.mat'

img = imgrr((86+15):(160+15),(96-15):(170-15),2:51);
img=mpi_upsampleImages(img);
img2 = img;
fimg = zeros(size(img));
[s1 s2] = size(img(:,:,1));
ww = fspecial('gaussian',[s1 s2],s1 /4);
varianceImage = zeros(s1,s2);
for i=1:s1
    for j=1:s2
        varianceImage(i,j) = var(img(i,j,:));
    end
end
ww = filter2(fspecial('gaussian',[s1 s2],s1/20),varianceImage.*ww/max(varianceImage(:)));
%ww = varianceImage.*ww./max(varianceImage(:));
figure(1);
surf(ww);
figure(2);
imagesc(ww);
figure(3);
surf(fspecial('gaussian',[s1 s2],s1 /8));


pause;

% simply refer each image to the last one
lastfImg = fft2(img(:,:,50).*ww);

%% use filter and shift things
displacements = zeros(2, 51);
Similarity = zeros(1,51);
for t=(51-2):-1:2
    current = fft2(img(:,:,t).*ww);
    auto = fftshift(ifft2(lastfImg.*conj(current)));
    figure(2)
    imagesc(img(:,:,t).*ww);
    axis image
    colormap gray
    figure(3)
    imagesc(abs(auto));
    axis image
    colormap gray
    pause(.1);
    Similarity(t) = max(auto(:));
    [row col] = find(auto == Similarity(t));
    displacements(:,t) = [row col] - [s1+1 s2+1]./2;
    magnitude = abs(displacements(:,t));
    if magnitude > 10
        disp(['We have a bad frame @ ' num2str(t)]);
    end
    img(:,:,t) = circshift(img(:,:,t),displacements(:,t));
    current = fft2(img(:,:,t).*ww);
    lastfImg = current;
end
figure(6);
plot(Similarity);

load '/v/raid1/npack/Processing/P043009_nufft/Output/endo_polyCoords.study27.slice1'
load '/v/raid1/npack/Processing/P043009_nufft/Output/epi_polyCoords.study27.slice1'


for t=2:(51-1)
    figure(7);
    subplot(1,2,1);
    imagesc(img2(:,:,t));
    hold on
    plot(endo_polyCoords_study27(:,1),endo_polyCoords_study27(:,2));
    plot(epi_polyCoords_study27(:,1),epi_polyCoords_study27(:,2));
    colormap gray
    subplot(1,2,2);
    imagesc(img(:,:,t));
    hold on
    plot(endo_polyCoords_study27(:,1),endo_polyCoords_study27(:,2));
    plot(epi_polyCoords_study27(:,1),epi_polyCoords_study27(:,2));
    colormap gray
    pause(.3);
end

return

n = 3;  % must be odd
circlist = zeros(size(img,1),size(img,2),n);
for i=1:n
    circlist(:,:,i) = fftshift(fft2(img(:,:,i)));
end
autoCorrelation = zeros([img(:,:,1) n n]);
DisplacementVectors = zeros([2 n n]);
for t = n:(68-n)
    fimg(:,:,t) = fftshift(fft2(img(:,:,t)));
    circlist(:,:,mod(t,n)+1) = fimg(:,:,t);
    for k=1:n
        i = k;
        j = n-k;
        autoCorrelation(:,:,i,j) = ifft2(circlist(:,:,mod(t,n)+1) .*  circlist(:,:,k));
        
        DisplacementVectors(:,i,j) = max(autoCorrelation(:,:,i,j));
    end
    
    figure(1);
    imagesc(img(:,:,t));
    axis image
    colormap gray
    figure(2)
    imagesc(log(abs(fimg(:,:,t))));
    axis image
    colormap gray
    pause(.1);
end
