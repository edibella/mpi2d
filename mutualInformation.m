function shifts_out = mutualInformation(uncroppedmovie, rangex, rangey,referenceFrame,sigma)
        
    if(exist('sigma') ~= 1)
        sigma = 1.27;
    end
    [sx sy tmax] = size(uncroppedmovie);
    
    %smooth the movie  Should really be doing anisotropic difussion or some other edge preserving filter
    %also we need it discritized to 256 so histogramming is simple
    gaussian = fspecial('gaussian',[sx sy],sigma);
    smoothed_movie = zeros(size(uncroppedmovie));
    for t=1:tmax
        smoothed_movie(:,:,t) = filter2(gaussian,uncroppedmovie(:,:,t));
    end
    refrenceImage = smoothed_movie(rangex,rangey,referenceFrame);
     
    %construct the different shift types
    granularityTypes = [1];
    limits = 10;
    possibleShifts = containers.Map({0},{[0 1 2 3 4]}); remove(possibleShifts,0);
    for i=1:length(granularityTypes)
        range = -limits:granularityTypes(i):limits;
        [X Y] = meshgrid(range,range);
        possibleShifts(granularityTypes(i)) = round([X(:) Y(:)]);
    end
    N = 256;
    h=zeros(N,N);
    shifts_out = zeros(tmax,2);
    % go through the movie, apply the shift and compute the NMI based on the 2d histogram
    %waith = waitbar(0,'fixing...');
    %finewait = waitbar(0,'checking all shifts');
    for t=1:(tmax-1)
        %waitbar(t/tmax,waith);
        disp(num2str(t));
        %figure(1);
        %differenceImage = abs(smoothed_movie(:,:,t) - smoothed_movie(:,:,t+1));
        differenceImage = abs(smoothed_movie(rangex, rangey,t) - refrenceImage);
        mymin = min(differenceImage(:)); mymax = max(differenceImage(:));
        %if(t == referenceFrame) mymin = 0; mymax = 1; end;
        %subplot(1,2,1), imagesc(differenceImage,[mymin mymax]), colormap gray;
        %title('before');
        for granularity = granularityTypes
            shifts = possibleShifts(granularity);
            clear NMI
            for shifti = 1:length(shifts)
                %waitbar(shifti/length(shifts),finewait);
                shift = shifts(shifti,:);
                shiftedImage = circshift(smoothed_movie(:,:,t),shifts(shifti,:));
                shiftedImage = shiftedImage(rangex, rangey);
                %compute the 2d histogram
                %h = jointhist(shiftedImage,refrenceImage,256);
                [h,C] = hist3(horzcat(shiftedImage(:),refrenceImage(:)),[256,256]);
%                     [r1 r2]=size(shiftedImage);
%                     Aratio = N/max(shiftedImage(:));
%                     Bratio = N/max(refrenceImage(:));
%                     tmat1 = ceil(shiftedImage*Aratio);
%                     tmat2 = ceil(refrenceImage*Bratio);
%                     h(:) = 0;
%                     for i=1:r1*r2
%                         h(tmat1(i),tmat2(i))=h(tmat1(i),tmat2(i))+1;
%                     end
                Pfr=h./sum(h(:));   %%%% Joint Image Intensity distribution.
                Pf=sum(Pfr,2); %%% Marginal Intensity distribution.
                Pr=sum(Pfr);
                temp = Pf(Pf~=0);
                Hf = -sum(temp.*log2(temp)); %marginal entropy for image 1
                temp = Pr(Pr~=0);
                Hr = -sum(temp.*log2(temp)); %marginal entropy for image 1
                temp = Pfr(Pfr~=0);
                Hfr = -sum(sum(temp.*(log2(temp)))); % joint entropy
                %compute NMI
                NMI(shifti)=(Hf+Hr')./Hfr;
            end
            [maxNMI, index] = max(NMI);
            shifts_out(t,:) = shifts_out(t,:) + shifts(index,:);
            smoothed_movie(:,:,t) = circshift(smoothed_movie(:,:,t),shifts(index,:));
        end
        %subplot(1,2,2), imagesc(abs(smoothed_movie(:,:,t) - refrenceImage),[mymin mymax]), colormap gray;
        %title('after');
        %pause(.1);
        
    end
    %close(waith);
    %close(finewait);
end