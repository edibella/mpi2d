
worked = -1;
whatDidNotWork = {};
parfid = fopen(ParFileName);
if parfid < 0
    whatDidNotWork(1) = {'Openning The File'};
    return;
end
parline = fgetl(parfid);
while ischar(parline)
    %disp(parline);
    %check for expressions that should evaluate to strings
    if ~isempty(strfind(parline,'infile')) || ~isempty(strfind(parline,'inendi')) || ~isempty(strfind(parline,'timeStampFile')) || ~isempty(strfind(parline,'intype'))
       tempIndex = strfind(parline,'=');
       parline = [parline(1:tempIndex) char(39) parline((tempIndex+1):end) char(39)];        
    end
    try
        eval([parline ';']);
    catch
        whatDidNotWork(end+1) = {parline};
    end
    parline = fgetl(parfid);
end
fclose(parfid);

worked = 1;