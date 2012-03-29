% param(1) = Cp0
% param(2) = Delta
% param(3) = alpha
% param(4) = tau
%
function f = Cpg(param,t)

Cp0 = abs(param(1));
Delta = param(2);
alpha = max(abs(param(3)),1);
tau = abs(param(4));

tprime = t-Delta;

alphaprime = alpha-1;
norm = exp(-alphaprime)*(tau*alphaprime)^alphaprime;

if (length(t) == 1)
    if (t < Delta) 
        f = 0;
    else
        f = (Cp0/norm)*(tprime^alphaprime)*exp(-tprime/tau);
    end;
else
    idx = find(t>=Delta);
    f = zeros(size(t));
    f(idx) = (Cp0/norm)*(tprime(idx).^alphaprime).*exp(-tprime(idx)/tau);
end;

return;