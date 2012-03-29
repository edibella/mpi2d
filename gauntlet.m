function stats = gauntlet()
load testbench.mat
for file=1:length(testfiles)
    ranget = rangets(file,1):rangets(file,2);
    auto = auto_registration_trial1(char(testfiles(file)),ranget,1);

    toexpress = ['true = offset.file' num2str(file) ';'];
    eval(toexpress);
    figure(file);
    plot(auto(:),true(:),'s');
    norms = sqrt(auto(:,1).*auto(:,1) + auto(:,2).*auto(:,2));
    axes
    
    %calculate r^2
    xmean = mean(auto(:));
    ymean = mean(true(:));
    sxx = sum((auto(:) - xmean).^2);
    syy = sum((true(:) - ymean).^2);
    sxy = sum((true(:) - ymean).*(auto(:) - xmean));
    stats.r2(file) = sxy/sqrt(sxx*syy);
    title(['Goodness = ' num2str(stats.r2(file))]);
    disp(['Goodness = ' num2str(stats.r2(file))]);
    temp = auto - true;
    stats.AverateError(file) = sum(sqrt(temp(:,1).^2 + temp(:,2).^2))/size(temp,1);
    disp(['Average |error| in pixels= ' num2str( stats.AverateError(file) )]);
    pause(.1);
end

end
   
