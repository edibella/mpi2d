function Ct = SBCmodelCtPkFixed(par,t)

global PKpar satThresh
 
Ncurves = size(PKpar,1);

Ct(1,:)=satBlood(par,t);

for ii=1:length(Ct(1,:))
    if(isnan(Ct(1,ii)))
        Ct(1,ii)=Ct(1,ii-1);
    end
end

% for i=2:Ncurves  %EVRD - needed??
%     tmp=sbc_ct(PKpar(i,:),[cpt(par,t),t]);
%     for ii=1:length(tmp)
%         if(isnan(tmp(ii)))
%             f=2;
%         end
%     end
% end

for i=2:Ncurves

    Ct(i,:) = sbc_ct(PKpar(i,:),[cpt(par,t),t]);

end;

 
return;

