%all_AstellasPatientsMain

%addpath /Users/ed/Documents/ShareWithWindows/matlab/mpi2d_Code_fromBrianAug2_2010

% scanner='Trio'
% % 'P043009','P061109',
%  for sPatient = {'P070909','P071709','P081009', 'P081309', 'P081809', 'P090309', 'P111309', 'P111609', 'P120109', 'P120409', 'P120709', 'P121509', 'P121609', 'P010710' }  % older ones
%    % tmpstring=char(strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',sPatient,'/Processing'));
%     cd(char(strcat('/v/raid1/ed/MRIdata/Cardiac/',scanner,'/',sPatient,'/Processing')));     delete('mpi2d.par')
%  parfiles = dir('*.par');
% 
% for p=1:length(parfiles)
%     if(strcmp(parfiles(p).name,'mpi2d.par')) continue; end    %kill off mpi2d
%     if(~isempty(strfind(parfiles(p).name,'rad'))) continue; end   %kill off subsets
%     
%     %extract the series and slice number
%     A = sscanf(parfiles(p).name,'series%d.slice%d.par')
%     if(isempty(A)) continue; end
%     series = A(1); slice = A(2); %subset = A(3);
% %     if(isempty(find(slicesToAllow == slice)))
% %         continue;
% %     end
%     
%     
%     copyfile(char(parfiles(p).name),'mpi2d.par');
%     mpi2d stg=[4.31]
% 
%     
% %     %find the delteSIcurves file
% %     template = strrep(strrep(parfiles(p).name,'series','study'),'par','mat');
% %     curvesCandidates = dir(['Output/deltaSIcurves*' template]);
% %     
% %         load(['Output/' curvesCandidates(1).name]);
%        
%        
% end
% end


scanner='Verio'
scanner='Trio'

%for sPatient = {'P051310','P051410','P052710','P062810'}
%     {'P050610', 'P040110','P041910',...
%     'P050510','P052010A',...
%     'P052010B','P060710','P061010'}
for sPatient = {'P071709'}

%     tmpstring=char(strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',sPatient,'/Processing'));
%     tmpstring_dest=char(strcat('/v/raid1/ed/MRIdata/Cardiac/',scanner,'/',sPatient));
    cd(char(strcat('/v/raid1/ed/MRIdata/Cardiac/',scanner,'/',sPatient,'/Processing')));     delete('mpi2d.par')
    
    parfiles = dir('*.par');
%curveFileNames = {};

for p=1:length(parfiles)
    if(strcmp(parfiles(p).name,'mpi2d.par')) continue; end    %kill off mpi2d
    if(~isempty(strfind(parfiles(p).name,'rad'))) continue; end   %kill off subsets
    
    %extract the series and slice number
    A = sscanf(parfiles(p).name,'series%d.slice%d.par')
    if(isempty(A)) continue; end
    series = A(1); slice = A(2); %subset = A(3);
%     if(isempty(find(slicesToAllow == slice)))
%         continue;
%     end
    
    
    copyfile(char(parfiles(p).name),'mpi2d.par');
    mpi2d stg=[4.31]

    
%     %find the delteSIcurves file
%     template = strrep(strrep(parfiles(p).name,'series','study'),'par','mat');
%     curvesCandidates = dir(['Output/deltaSIcurves*' template]);
%     
%         load(['Output/' curvesCandidates(1).name]);
       
       
end


    
end




% looking at one by one,  8/5/10  EVRD, to select which slice to pull AIF
% from. 
% All multi-srt: 
% verio:
% P040110  looks like 2 slices per injection. rest: use 18000 (slice 1)
% adeno: slice 1 (22000)
% lexi slice 1 (67000)
% 
% 
% P041910:
% wow, 4 slices each,
% rest: 22000 to 25000 1-4 pretty close except slice 3 low. 
% adeno: like slice 2, though #1 is highest
% lexi: looks like problem with 74-77,
% 79000 (slice 2) pretty good, again #1 workable.. 
% 
% P050510:
% 3 slice each,
% rest: really differnt, slice 1 highest. 
% adeno: again very different, slice 1 
% lexi: not so different, slice 1
% 
% P052010A:
% 3 slices each,
% rest: slice 1 way higher, 2-3 similar
% adeno: 1 highest
% lexi: 1 high
% 
% P052010B: 
% all 3 slices but lexi only 2. 
% slice 1 for all
% 
% P060710:
% 4 slices rest, 3 adeno and lexi
% slice 1 for all.. 
% 
% P061010:
% 3 rest, 2 adeno and lexi.
% All use slice 1. Amazing how much lower slice 2 is! (2 and 3 for rest)
% hmm, looked with osirix and seemed lexi 64000 and 65000 AIF looks similar on images... 
%     the tails do match, but peak almost double... will check with 3.11