function shifts = GaneshModel(fullcinemri,rangex,rangey,ranget,RVpoint,referenceFrame)
if(~exist('RVpoint'))
    [RVpoint, ~] = findRVLV(cinemri1);
end
oldpwd = pwd;
addpath('/v/raid1/bmatthew/registration_code');
cd('/v/raid1/bmatthew/registration_code');
%programatically perform stage 3.1 by creating the necessary files for
%stage 3.2 to run
RVbloodContour = zeros(12,2);
index = 1;
for theta = 0:(2*pi/12):(2*pi)
    RVbloodContour(index,1) = 3*sin(theta) + RVpoint(1);
    RVbloodContour(index,2) = 3*cos(theta) + RVpoint(2);
    index = index + 1;
end
epiContour = zeros(12,2);
index = 1;
for theta = 0:(2*pi/12):(2*pi)
    epiContour(index,1) = 5*sin(theta) + 20;
    epiContour(index,2) = 5*cos(theta) + 20;
    index = index + 1;
end
endoContour = zeros(12,2);
index = 1;
for theta = 0:(2*pi/12):(2*pi)
    endoContour(index,1) = 3*sin(theta) + 20;
    endoContour(index,2) = 3*cos(theta) + 20;
    index = index + 1;
end
angle = 0; %#ok<NASGU>
save('Output/blood_polyCoords.study12.slice3','RVbloodContour','-ascii');
save('Output/epi_polyCoords.study12.slice3','epiContour','-ascii');
save('Output/endo_polyCoords.study12.slice3','endoContour','-ascii');
save('Output/Roi_start_angle.study12.slice3','angle','-ascii');
%save the cine file somewhere  probably where the infile is pointing
if(length(ranget) ~= size(fullcinemri,3))
    cinemri = fullcinemri(rangex,rangey,ranget);
else
    cinemri = fullcinemri(rangex,rangey,:);
end
cinemri1 = cinemri;
save('Output/cinemri.study12.slice3.mat','cinemri');
save('Output/cinemri1.study12.slice3.mat','cinemri1');
if(length(ranget) ~= size(fullcinemri,3))
    buffer_cinemri = fullcinemri(:,:,ranget);
else
    buffer_cinemri = fullcinemri(:,:,:);
end
save('Output/buffer_cinemri.mat','buffer_cinemri');
clear temp

%temp.rangex = max(min(rangex),1):min(max(rangex),size(fullcinemri,1));
%temp.rangey = max(min(rangey),1):min(max(rangey),size(fullcinemri,2));
temp.rangex = [num2str(min(rangex)) ':' num2str(max(rangex))];
temp.rangey = [num2str(min(rangey)) ':' num2str(max(rangey))];
temp.ranget = [num2str(min(ranget)) ':' num2str(max(ranget))];
temp.referenceFrame = referenceFrame;
updateParFile('mpi2d.par',temp);



%here's ganesh's code
mpi_model_reg;


%here's the shifts he leaves behind
myshifts = load('g_sh.mat');
rmpath('/v/raid1/bmatthew/registration_code');
cd(oldpwd);
shifts = -squeeze(myshifts.g_sh(:,:,4));