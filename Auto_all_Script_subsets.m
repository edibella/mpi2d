%make sure that all mother data sets have been processed to the point of
%contour creation at least.
mpi2d stg=0.3   %for each subset copy shifts, timestamp and contour data from mother to child
if(exist('mpi2d.par')) delete('mpi2d.par'); end
files = dir([folder '/rad.*.par']);
for i=1:length(files)
    copyfile(char(files(i).name),'mpi2d.par');
    disp(['i = ' num2str(i)]);
    disp(['Processing ' files(i).name]);
    mpi2d stg=[0.23 3.2 3.4]
end
return;