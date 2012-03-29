% correct built in gammainc for negative x where abs(x) > a+1
function out = gammaincc(x,a)

out = zeros(size(x));

f = find(x > -(a+1));

out(f) = gammainc(x(f),a);

f = find(x <= -(a+1));

for i=1:length(f)
    idx = f(i);
%     out(idx) = gammainc(x(idx),a);
    out(idx) = gammastarseries(x(idx),a);
end;

return;

function gstar = gammastarseries(x,a)

% checked vs. Mathematica Gamma[a,0,x]
% Abramowitz and Stegun Eq. 6.5.4 and 6.5.29
gstar = 1;
sum = 0;
n = 0;

while (abs(gstar) > 1e-10 & n <= 4*abs(x))
    tmp1 = n*log(x);
    tmp2 = exp(tmp1-x);
    tmp3 = gamma(a+n+1);
    gstar = tmp2/tmp3;
    sum = sum+gstar;
    n=n+1;
end;

gstar = real(sum*x^a);

% % use recurrence - bug somewhere
% % n = 0 term
% gstar = exp(-x)/gamma(a+1);
% sum = gstar;
% n = 1;
% 
% while (abs(gstar) > 1e-10 & n <= 4*abs(x))
%    gstar = gstar*x/(a+n+1);
%    sum = sum + gstar;
%    n = n+1;
% end;
% 
% g = real(sum*x^a);

return;
