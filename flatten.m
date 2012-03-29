function ret = flatten(in)

ret = reshape(in,[1 prod(size(in))]);

return;