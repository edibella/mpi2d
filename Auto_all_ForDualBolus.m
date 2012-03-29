%1. create a folder to hold the patient data and processed output data
%2. cd to that folder in matlab
%3. pop in the CD from the scanner, or if you already have the dicom data in a
%   folder skip this step.  The dicoms will be retrieved from there in step 4
%4. Run mpi2d stg=1.1 and point it to the CD folder that contains the dicoms
%       at this point you may either delete folders you don't want to process, or
%       you can do this filtering later
%5. Run mpi2d stg=1.2
%6. If you want to omit the processing of certain slices then delete the
%       corresponding par files now.  They are located in the processing directory(current)
%7. Change the folder path in the line below to point to the Processing folder that it
%       should already be in.  (his should match the current directory but in case it does not, 
%       setting the folder this way ensures that this dataset will be processed)
folder = pwd;
%folder = '/home/mirl/bmatthew/Desktop/test/Files you should get/Processing';
%8. Run this script

first_series=52000

cd(folder);

for islice=4:5
    iseries=first_series+(islice-1)*1000
    for isub=90:93  % low dose sys, low dose dias, high sys, high dias
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
      mpi2d stg=[2.8 3.12 3.2 3.4 4.1] myothreshold=.5
    end
end
return;

%------------------------------
%the following snippets show the recommended order of fixing datasets 

%for datasets that didn't segment well, tell mpi2d which dataset to do
copyfile('seriesX,sliceY.par','mpi2d.par');  %replace X and Y with those of the dataset you wish to modify

%------------------------------
%then look at the weighted myocardium image (typically on the second figure for that dataset
%then choose a different threshold and run
mpi2d stg=[3.12 3.2 3.4 4.1 6.1] myothreshold=.6

%------------------------------
%if it still has a hard time or you would just like to do it yourself
%un-comment the following line if you want to zero out the shift vector ie. the automatic registration completely botched it
%mpi2d stg=2.11
mpi2d stg=[3.11]
%then once you're done
mpi2d stg=[0.23 3.2 3.4 4.1 6.1]

%------------------------------
% This section shows an example of how to tweak the final k-trans display
% stage.  The Alpha is how much of the myocardium's k-trans color should
% show up, and the MyoFatness is how wide to make the tappered ring
mpi2d stg=6.1 Alpha=.2 MyoFatness=.1

%------------------------------
%This following section is for complete manual registration and
%segmentation of all datasets.

cd(folder);
if(exist('mpi2d.par')) delete('mpi2d.par'); end
files = dir([folder '/*.par']);
for i=1:length(files)
    copyfile(char(files(i).name),'mpi2d.par');
    disp(num2str(i));
    disp(files(i).name);
    %un-comment the following line if you want to zero out the shift vector ie. the automatic registration completely botched it
    %mpi2d stg=2.11
    mpi2d stg=[3.11]
    keyinput = input('press 1');
    while(keyinput ~= 1)
        pause(1);
    end
    %0.23 applied the shifts that exist so that there is no image wrapping
    %the other stages are described above
    mpi2d stg=[0.23 3.2 3.4 4.1 6.1]
end

