function addImageToGauntlet(filename)
load('testbench.mat');
testfiles(end+1) = {filename};
load(filename);
disp('please outline a cropping');
figure(10),imagesc(cinemri(:,:,5));colormap gray;
[BW x y] = roipoly;
if(~isempty(x))
    for t=1:size(cinemri,3)
        temp(:,:,t) = cinemri(min(x):max(x),min(y):max(y),t);
    end
    cinemri = temp;
end
save(filename,'cinemri');
rangets = vertcat(rangets,[2 size(cinemri,3)]);
ranget = 2:size(cinemri,3);
true = auto_registration_trial1(filename,ranget,1);
toexpress = ['offset.file' num2str(size(testfiles,2)) ' = true;'];
eval(toexpress);
save('testbench.mat','offset','testfiles','rangets');
end