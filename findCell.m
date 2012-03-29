function index = findCell(array,value)
index = find(cellfun(@(y) any(strcmp(y,value)),array));
end