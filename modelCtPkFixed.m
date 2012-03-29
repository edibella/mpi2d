function Ct = modelCtPkFixed(par,t)

global PKpar
 
Ncurves = size(PKpar,1);

 
for i=1:Ncurves

    Ct(i,:) = sbc_ct(PKpar(i,:),[cpt(par,t),t]);

end;

 
return;

