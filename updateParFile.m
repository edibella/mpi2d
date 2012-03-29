function updateParFile(filename, UpdatingstructWithStringFields)
thingsToChange = fieldnames(UpdatingstructWithStringFields);
fid = fopen(filename);
outputlines = {};
tline = fgetl(fid);
while ischar(tline)
    trimmed = strtrim(tline(1:(strfind(tline,'=')-1)));
    for i=1:length(thingsToChange)
        if(strcmp(trimmed,thingsToChange{i}))
            if(~ischar(UpdatingstructWithStringFields.(thingsToChange{i})))
                tline = [thingsToChange{i} ' = ' num2str(UpdatingstructWithStringFields.(thingsToChange{i})) ';'];
            else
                tline = [thingsToChange{i} ' = ' UpdatingstructWithStringFields.(thingsToChange{i}) ';'];
            end
            break;
        end
    end
    outputlines{end+1} = tline;
    tline = fgetl(fid);
end
fclose(fid);

fid = fopen(filename,'w');
for i=1:length(outputlines)
    fprintf(fid,'%s\n',outputlines{i});
end
fclose(fid);