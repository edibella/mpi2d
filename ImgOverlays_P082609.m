

% %%%%% FOR PIXELWISE DATA
% 
% % First to loadthe panel images for this pt....
% load /v/raid1/npack/Processing/P082609_nufft/Output/PanelImages_P082609_ser33.mat
% Img=PanelImages_P082609_ser33(:,:,16); % this reference image was manually selected
% % Next to load the Ktrans overlays for this pt....
% load A.mat; load B.mat; load C.mat
% KtransMask=[C B A];
% 
% %%% To create a pixelwise transparency map for overlaying the perfusion mask
% TransParencyMask=zeros(size(KtransMask));
% TransParencyMask(find(KtransMask))=0.5;  % This const. value can be changed to alter the transparency of the overlay image
% 
% spect=load('spect.cmap');
% m=max(max(spect));
% spect=spect/m;
% 
% 
% figure
% subimage(1) % so that the cardiac image has it's own colormap
% A=imagesc(Img,[0 max(Img(:))]);
% axis image,hold on,colormap gray
% subimage(2) % so that the overlay image has it's own colormap
% H=imagesc(KtransMask,[0 1]);colormap(spect)
% set(H,'AlphaData',TransParencyMask),colorbar










%%%%% FOR PIXELWISE DATA

% First to load the panel images and the Ktrans overlays for this pt....
load /v/raid1/npack/Processing/P082609_nufft/Output/PanelImages_P082609_ser33.mat
Img=PanelImages_P082609_ser33(:,:,16); % this reference image was manually selected
load /v/raid1/npack/Processing/P082609_nufft/A.mat; 
load /v/raid1/npack/Processing/P082609_nufft/B.mat; 
load /v/raid1/npack/Processing/P082609_nufft/C.mat
KtransMask=[C B A];



figure
subimage(1) % so that the cardiac image has it's own colormap
A=imagesc(Img,[0 max(Img(:))]);  axis image,hold on,colormap gray


%%% embed colorbar...-----------------------------------------------------
subimage(2)
maxx=max(KtransMask(:));
minn=0.8*min(KtransMask(:)); %minn=0   % can comment this if want range min to max

nrows=size(KtransMask,1)
delta_signal=(maxx-minn)/nrows

vec=(nrows:-1:1)* delta_signal + minn;
KtransMask(:,end-5:end)=[vec; vec; vec; vec; vec; vec;]';  % replace far right columns with "colorbar"
%%% To create a pixelwise transparency map for overlaying the perfusion mask
TransParencyMask=zeros(size(KtransMask));
TransParencyMask(find(KtransMask))=0.5;  % This const. value can be changed to alter the transparency of the overlay image
spect=load('spect.cmap');  m=max(max(spect));  spect=spect/m;
H=imagesc(KtransMask,[0 1]);colormap(spect)

% label it 
pp=gca;   set(pp,'YAxisLocation','right')

% choose number of tick marks and their location, based on row number (ie. try 11 evenly spaced)
YlimRange=get(pp,'Ylim');
YlimRange=YlimRange(2)-YlimRange(1);
tickLocs=(0:10)*(YlimRange/10)+0.5;
set(pp,'YTick',tickLocs);
tickLabels=(10:-1:0)*((maxx-minn)/10)+ minn;    % not sure will line up exactly, may need to tweak  (e.g. change 0.5 above)
set(pp,'YTickLabel',tickLabels)
%set(pp,'YTickLabel',[5 4 3 2 1 0 -1 -2 -3 -4 -5])
%%%%----------------------------------------------------------------------

set(H,'AlphaData',TransParencyMask)%,colorbar














% 
% 
% %%%%% FOR PIXELWISE DATA
% 
% % First to loadthe panel images for this pt....
% load /v/raid1/npack/Processing/P082609_nufft/Output/PanelImages_P082609_ser33.mat
% Img=PanelImages_P082609_ser33(:,:,16); % this reference image was manually selected
% % Next to load the Ktrans overlays for this pt....
% load Apixel.mat; load Bpixel.mat; load Cpixel.mat
% KtransMask=[Cpixel Bpixel Apixel];
% 
% %%% To create a pixelwise transparency map for overlaying the perfusion mask
% TransParencyMask=zeros(size(KtransMask));
% TransParencyMask(find(KtransMask))=0.5;  % This const. value can be changed to alter the transparency of the overlay image
% 
% spect=load('spect.cmap');
% m=max(max(spect));
% spect=spect/m;
% 
% 
% figure
% subimage(1) % so that the cardiac image has it's own colormap
% A=imagesc(Img,[0 max(Img(:))]);
% axis image,hold on,colormap gray
% subimage(2) % so that the overlay image has it's own colormap
% H=imagesc(KtransMask,[0 1]);colormap(spect)
% 
% 
% 
% 
% set(H,'AlphaData',TransParencyMask),colorbar