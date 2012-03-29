% param(1) = Cp0
% param(2) = delta
% param(3) = alpha
% param(4) = tau
% param(5) = T
%
function f = Cps(param,t)

Cp0 = abs(param(1));
Delta = param(2);
alpha = max(abs(param(3)),1);
tau = abs(param(4));
T = abs(param(5));

if(~isreal(Delta))
    Delta=abs(Delta);
end
tprime = t-Delta;

norm = gamma(alpha);

xi = 1/tau-1/T;
gam = norm;

if alpha<1
    alpha=1;
end

if (length(t) == 1)
    if (t < Delta) 
        f = 0;
    else
        f = (Cp0/norm)*gam*exp(-tprime/T)*gammaincc(xi*tprime,alpha);
    end;
else
    idx = find(t>=Delta);
    f = zeros(size(t));
%     if(~isreal(xi*tprime(idx)) | ~isreal(alpha))
%         disp(['xi: ', num2str(xi), ' alpha: ',num2str(alpha), ' delta: ',num2str(Delta)])
%         Delta=abs(Delta);
%         disp(['Delta: ',num2str(Delta)]);
%         pause
%     end
    f(idx) = (Cp0/norm)*gam*exp(-tprime(idx)/T).*gammaincc(xi*tprime(idx),alpha);
end;

return;