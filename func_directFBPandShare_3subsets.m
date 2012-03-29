function img=func_directFBPandShare_3subsets(kSpace, template_dicom_directory, output_dicom_directory, slicenum, seriesnum, seriesdesc )

% Program to reconstruct images from radially sampled MRI data using Complex Filtered Backprojection

% Revisiting 7/16/07, EVRD.  Old version was to just recon a single frame,
% with 96 projections split into 4 subsets, into 4 images. used STCR
% (wonder now if temporal sort of wrong due to slight rotations. If use
% correct angle each time, reconstructed images should not rotate. Not sure
% if they do rotate (wrong angle order?), or if just streaking rotates...
% order sort of looks like 1,3,2,4 based on SRT...

% Now wish to use with new data where get offset rotation each time frame (period 
% of 4, but order is 1,3,2,4. Order within a time frame is also 1,3,2,4)
% so now really can do STCR (and sliding window, etc.)  
% This new data is with different SRTs in different slices...

% revisiting again, EVRD 12/28/09 or so. Want to use for 72 rays, 3 subsets
% (though acquired with 12 subsets of 6 rays each.).   For now, don't use
% neighbors in time (see how that is done in Eugene's code, which is
% heavily borrowed here. 




%% Flags

% InterleaveSampling=1 - projections for adjacent time frame are interleaved
InterleaveRecon = 0;
InterleaveSampling = 1;
interleaveFactor = 4;
if InterleaveSampling == 0
   interleaveFactor = 1;
   InterleaveRecon = 0;
end
InterleaveSampling = 1 %1;
% SubsetRecon = 1 - only projection of the current subset included in recon
SubsetRecon = 0 %1;

% LowResRecon = 1 - high frequency reduced to achieve streaking-less recon
LowResRecon = 0;

Nr=256


nn = ndims(kSpace)

%nRead, nAngles, nCoils, nFrames, nSlices 
[Nx, Ny, Nc,Nt, Ns] = size(kSpace);
% First and Last Time Frames for Recon
Nts = 1;  %10;  % assumes first frame is 1
Ntf = Nt;  %10  %11;

nCoils=Nc; 
%   kSpace = single(reshape(kSpace,Nx,Ny,Nc,Nt,Ns));
   
interleavesFrameToFrame=4
interleavesWithinFrame=12
permuteIndex=[1 3 2 4];
permuteIndex=[1 2 3 4];

thetaMatrix = zeros(Ny,interleavesFrameToFrame);   %,'single');
for a=1:interleavesFrameToFrame
        theta = 0:1:Ny-1;
        theta = 180/Ny*theta + (a-1)*180/Ny/interleavesFrameToFrame;
        thetaMatrix(:,a) = reshape(theta,Ny,1);
end

%%
%% Partial Recon
% Parameter fraction defines how many projections will give low frequency contributions
% fraction can be equal to 0.25, 0.5, 1.0 
%fraction = 0.25
scaleCentralRays=3   % what to scale up by for central rays. Equal to 1/fraction
% Parameter offset defines projections from what subset give low frequency contributions

modifyWindowSize=0.8
   Nxc = Nx/2+1;
   if SubsetRecon == 0
      Nwf = round(65*modifyWindowSize);
      Nones = round(17*modifyWindowSize);
      N2 = (Nwf-Nones)/2;
      hmw0 = hanning(Nwf-Nones);
      hmw = ones(Nwf,1);
      hmw(1:N2) = hmw0(1:N2);
      hmw(end-N2+1:end) = hmw0(N2+1:end);
      wfHR = ones(Nx,1);
      wfHR(Nxc-(Nwf-1)/2:Nxc+(Nwf-1)/2) = wfHR(Nxc-(Nwf-1)/2:Nxc+(Nwf-1)/2)-hmw;
      wfLR = ones(Nx,1);
      wfLR(Nxc-(Nwf-1)/2:Nxc+(Nwf-1)/2) = wfLR(Nxc-(Nwf-1)/2:Nxc+(Nwf-1)/2)+hmw*(scaleCentralRays-1);
   end
   if SubsetRecon == 1 && LowResRecon == 1
      Nwf = Ny+1;
      Nones = 13;
      N2 = (Nwf-Nones)/2;
      hmw0 = hanning(Nwf-Nones);
      hmw = ones(Nwf,1);
      hmw(1:N2) = hmw0(1:N2);
      hmw(end-N2+1:end) = hmw0(N2+1:end);
      wfHR = zeros(Nx,1);
      wfLR = zeros(Nx,1);
      wfLR(Nxc-(Nwf-1)/2:Nxc+(Nwf-1)/2) = hmw;
   end
   if SubsetRecon == 1 && LowResRecon == 0
      wfLR = ones(Nx,1);
   end
      
figure(22); clf; plot(wfLR,'r')
hold on
plot(wfHR,'k')
legend('wfLR','wfHR')
disp('checking low and high weighting')
%pause

%s=1;
clear theta
nRaysPerSubset=Ny/interleavesWithinFrame;
nSubsets=12;
nSlices=Ns
for subset=1:3
for islice=1:1   %slices   %1:nSlices
    for icoil=1:nCoils
        for t=Nts:Ntf
            frameOrder = permuteIndex(mod(t-1,interleavesFrameToFrame)+1);  % goes 1,2,3,4
            index=0;
       %for offset=[2 8 5 11 1 7 4 10]; %0]  % 1 2 3]  %[0 2 1 3]   % same as 1 3 2 4 
           
       
            if subset==1
                for offset=[0 6 3 9]
                    index=index+1;
                    subsets_kSpace(:,:,index) = single(squeeze(kSpace(:,1+offset:interleavesWithinFrame:end,icoil,t,islice))).*repmat(wfLR,[1,nRaysPerSubset]);
;
                    theta(:,index)=squeeze(thetaMatrix(offset+1:interleavesWithinFrame:end,frameOrder));
                end
                for offset=[2 8 5 11 1 7 4 10]
                    index=index+1;
                    subsets_kSpace(:,:,index) = single( squeeze(kSpace(:,1+offset:interleavesWithinFrame:end,icoil,t,islice)) ).*repmat(wfHR,[1,nRaysPerSubset]);
;
                    theta(:,index)=squeeze(thetaMatrix(offset+1:interleavesWithinFrame:end,frameOrder));
                end
            elseif subset==2
                for offset=[2 8 5 11]
                    index=index+1;
                    subsets_kSpace(:,:,index) = single(squeeze(kSpace(:,1+offset:interleavesWithinFrame:end,icoil,t,islice))).*repmat(wfLR,[1,nRaysPerSubset]);
;
                    theta(:,index)=squeeze(thetaMatrix(offset+1:interleavesWithinFrame:end,frameOrder));
                end
                for offset=[0 6 3 9 1 7 4 10]
                    index=index+1;
                    subsets_kSpace(:,:,index) = single( squeeze(kSpace(:,1+offset:interleavesWithinFrame:end,icoil,t,islice)) ).*repmat(wfHR,[1,nRaysPerSubset]);
;
                    theta(:,index)=squeeze(thetaMatrix(offset+1:interleavesWithinFrame:end,frameOrder));
                end
            elseif subset==3
                for offset=[1 7 4 10 ]  %[0 6 3 9]
                    index=index+1;
                    subsets_kSpace(:,:,index) = single(squeeze(kSpace(:,1+offset:interleavesWithinFrame:end,icoil,t,islice))).*repmat(wfLR,[1,nRaysPerSubset]);
;
                    theta(:,index)=squeeze(thetaMatrix(offset+1:interleavesWithinFrame:end,frameOrder));
                end
                for offset=[0 6 3 9 2 8 5 11]
                    index=index+1;
                    subsets_kSpace(:,:,index) = single( squeeze(kSpace(:,1+offset:interleavesWithinFrame:end,icoil,t,islice)) ).*repmat(wfHR,[1,nRaysPerSubset]);
;
                    theta(:,index)=squeeze(thetaMatrix(offset+1:interleavesWithinFrame:end,frameOrder));
                end
            end
            
        
       theta_all(:,t)=reshape(theta,nSubsets*nRaysPerSubset,1);
       subsets_kSpace_all(:,:,icoil,t,islice)=reshape(subsets_kSpace,Nx,nSubsets*nRaysPerSubset,1);
        end

    end
end
   
%%
% REad in one time frame, and then have 4 subsets. Recon one of them, no
% constraints. Then recon all 4 together. 
   

for islice=1:1    %slices    %1:nSlices
    islice
  for icoil=1:nCoils
      icoil
  %  for i=4:nTimes
    for i=Nts:Ntf
        theta1=squeeze(theta_all(:,i));  % theta the same for all 
   %     tmp_kSpace=squeeze(kSpace(:,:,i,islice));
        img_est = single(iradonKC_Single( subsets_kSpace_all(:,:,icoil,i,islice),theta1,'linear','Ram-Lak',1,Nr));
%         figure(1),imagesc(abs(img_est(:,:))),colormap gray,brighten(0.3); title('direct recon')
%         title(strcat('Slice  ',num2str(islice),'  coil  ',num2str(icoil), '  time ',num2str(i)))
       % pause
        imgThisSlice(:,:,i,icoil) = rot90(img_est);  % rot90 added for P121609  EVRD 12/21/09
    end
  end

imgSqSofS(:,:,:,islice) = squeeze(sum(imgThisSlice.*conj(imgThisSlice),4));  % rot90 added for P121609  EVRD 12/21/09
   

 end
   
   imgSqSofS=sqrt(imgSqSofS) ;
   
   
   
   outfile=strcat(seriesdesc,int2str(subset),'_window',int2str(10*modifyWindowSize));
   
   mkdir(strcat(output_dicom_directory,int2str(seriesnum+subset),')/'))
   
outname = strcat(output_dicom_directory,int2str(seriesnum+subset),')/',outfile,'_');
%mkdir(strcat(output_dicom_directory,'/Recons_subset',int2str(subset)) )

       
  % outname=strcat(outname,int2str(slice),'_subset',num2str(subset),'_window',int2str(10*modifyWindowSize))
  % save(outfile, 'imgSqSofS')
  [sx sy sz] = size(imgSqSofS);
  cropsx=(sx/4+1):(sx-sx/4);
  cropsy=(sy/4+1):(sy-sy/4);
  cropsx=63:(63+128);
  cropsy=63:(63+128);
  %crops=1:Nr;
  
  for bb=1:size(imgSqSofS,3)  % start at frame 2 so matches the recons TaeHo and Nate make with TCR
   %imgrr(:,:,bb)=imgSqSofS(crops,crops,bb);   % CHECK SOMETIME IF OFF BY 1!!  %now putting back first frame...
   imgrr(:,:,bb)=flipud(imgSqSofS(cropsx,cropsy,bb));   % CHECK SOMETIME IF OFF BY 1!!  %now putting back first frame...
   %imgrr(:,:,bb-1)=5000*flipud(imgSqSofS(crops,crops,bb));   % CHECK SOMETIME IF OFF BY 1!!
  end
  
  seriesdesc_out = strcat(seriesdesc,int2str(subset),'_window',int2str(10*modifyWindowSize));
  
  % write out as dicom:
  % first skip the first frame to match TCR:
  imgrr=imgrr(:,:,2:end);  %EVRD 5/19/10
   %func_dicomHeader_noNormalizing(template_dicom_directory, imgrr, outname, seriesnum+subset, seriesdesc_out);
   func_dicomHeader_OrderByTime_SkipFirst(template_dicom_directory, imgrr, outname, seriesnum+subset, seriesdesc_out);
   % what happens if write out dicom with same series number, different
   % description and filename?   Or, how decide on series number, say 1st
   % subset is always 1000+series, second is 2000+series, and 3rd is 3000+
   % series
   clear imgrr
   
   

clear subsets_kSpace, subsets_kSpace_all;
end  % subset loop

  
