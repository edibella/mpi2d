function smoothed = anisotropicDifuse(img,k)
    smoothed = img;
    fil = zeros(3,3,8);
    fil(:,:,1) = [0 0 0;0 -1 1;0 0 0];
    fil(:,:,2) = fliplr(squeeze(fil(:,:,1)));
    fil(:,:,3) = flipud(squeeze(fil(:,:,1))');
    fil(:,:,4) = flipud(squeeze(fil(:,:,3)));
    fil(:,:,5) = [0 0 1;0 -1 0;0 0 0];
    fil(:,:,6) = fliplr(squeeze(fil(:,:,5)));
    fil(:,:,7) = flipud(squeeze(fil(:,:,6)));
    fil(:,:,8) = fliplr(squeeze(fil(:,:,7)));
    for t=1:20
        timg = zeros(size(img));
        for i=1:8
            f = filter2(fil(:,:,i),smoothed);
            c = exp(-(f/k).^2);
            timg = timg + 1/7*f.*c;
        end
        smoothed = smoothed + timg;
        figure(2)
        imagesc(smoothed);
        colormap gray;
        pause();
    end
end