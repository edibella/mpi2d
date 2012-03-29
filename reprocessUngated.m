function reprocessUngated(raid1path, first_series, numSlices, startSlice, preBolusFlag);
% to Process ungated perfusion.  9/19/11
% SET below parameters
% AFTER it is done, run 3.11 for all and fix, then mpi2d stg=[0.23 3.2 3.4 4.1
% 4.2]

% NOTE: Two functions put below in this file. Used to be separate scripts.


% reads in dicom like in 1.2, or
% read in cinemri, just recon with no header?
% find rv and lv centers.. 
% run self-gating

if preBolusFlag==1
    isub_range=[90 91]
else
    isub_range=[92 93]
end




for islice=startSlice:numSlices  % chnage also dir name below!!
iseries=first_series+(islice-1)*1000;

cd([raid1path,'/Processing/']);

% next part is from Auto_all_ForDualBolus
    
    for isub=isub_range  % low dose sys, low dose dias, high sys, high dias
      %  if isub==90 || isub==92
      parfile=['series',int2str(iseries+isub),'.slice',int2str(islice),'.par'] 
      if(exist('mpi2d.par')) delete('mpi2d.par'); end
      copyfile(parfile,'mpi2d.par');
    
      disp(['Processing ' parfile]);
    %2.6 cross-correlation registration with gaussian centering
        %if this registration screws up, run mpi2d stg=2.11 to apply zero
        %shifts then do 3.11 to do manual registration
    %3.12 automatic segmentation
    %3.2 extract curves
    %3.4 convert curves to deltaSIcurves
    %4.1 compute 2-compartment model params
    %6.1 display k-trans
      %mpi2d stg=[2.8 3.12 3.2 3.4 4.1] myothreshold=.5
      mpi2d stg=[0.23 3.2 3.4 4.1]
    end


end





