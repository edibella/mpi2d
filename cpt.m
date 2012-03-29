function Cp = cpt(par,t,fplot)

if nargin < 3
    fplot=0;
end

delta1 = par(1);
alpha1 = par(2);
tau1 = par(3);
A2 = par(4);
delta2 = par(5);
tau2 = par(6);
A3 = par(7);
delta3 = par(8);
A4 = par(9); 
T = par(10);
% if length(par) < 11
%     A1 = 1;
% else
    A1 = par(11);
%     A2 = A2*6;
%     A3 = A3*6;
%     A4 = A4*6;
% end
 
Cp = Cpg([A1 delta1        alpha1 tau1],t) + ...
     Cpg([A2 delta1+delta2 alpha1 tau2],t) + ...
     Cpg([A3 delta1+delta3 alpha1 tau2],t) + ...
     Cps([A4 delta1+delta3 alpha1 tau2 T],t);


% if length(par) ==10
%     Cp = 6.0*Cp;
% end

if(fplot)
    plot(t,Cpg([A1 delta1        alpha1 tau1],t),'r-.')
    hold on
    plot(t,Cpg([A2 delta1+delta2 alpha1 tau2],t),'b-.')
    plot(t,Cpg([A3 delta1+delta3 alpha1 tau2],t),'g-.')
    plot(t,Cps([A4 delta1+delta3 alpha1 tau2 T],t),'c-.')
    plot(t,Cp,'k-','LineWidth',2)
    hold off
end
    
 
% Cp = 6.0*(mean(exp(-(t(last)-delta1)/14.0))/mean(Cp(last)))*Cp;
% Cp = (0.5/Cp(length(t)))*Cp;
Cp=real(Cp);
 
return;

