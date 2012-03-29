folder = '/v/raid1/bmatthew/MRIdata/Cardiac/Verio/P041910/Processing';
cd(folder);
delete('mpi2d.par');
files = dir([folder '/*.par']);
% for i=55:57
%     copyfile(char(files(i).name),'mpi2d.par');
%     disp(files(i).name);
%     mpi2d stg=[3.11]
%     keyinput = input('press 1');
%     while(keyinput ~= 1)
%         pause(1);
%     end
% end

for i=55:57
    copyfile(char(files(i).name),'mpi2d.par');
    mpi2d stg=[3.2 3.4 4.1]

end
