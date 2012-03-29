folder = '/v/raid1/bmatthew/MRIdata/Cardiac/Verio/P060710/Processing';
cd(folder);
if(exist('mpi2d.par')) delete('mpi2d.par'); end;
files = dir([folder '/*.par']);
for i=12:length(files)
    copyfile(char(files(i).name),'mpi2d.par');
    disp(files(i).name);
    ParFileName = files(i).name;
    mpi2d stg=[3.11]
    keyinput = input('press 1');
    while(keyinput ~= 1)
        pause(1);
    end
    mpi2d stg=0.23
    mpi2d stg=[3.2 3.3 3.4 4.1 4.2]
end