% Tissue function

function f=sbc_ct(kpar,aift)


if (size(aift,1) > size(aift,2)); aift = aift'; end;
    
N=length(aift);

aif=aift(1:N/2);
t=aift(N/2+1:end);
lent=length(t);

delta_t=(t(lent)-t(1))/lent;
sh_aif=interp1(t,aif,t-kpar(4),'cubic');


% f=ifft(fft(kpar(1)*aif).*fft(exp(-kpar(2)*t)));
% y=kpar(1)*exp(-kpar(2)*t);
% 
% X=fft([aif zeros(1,length(y)-1)]);
% Y=fft([y zeros(1,length(aif)-1)]);
% 
% f=ifft(X.*Y);

% and now with the convolution directly

% f=conv(aif,kpar(1)*exp(-kpar(2)*t));
f=conv(sh_aif,delta_t*kpar(1)*exp(-kpar(2)*t));

f=f(1:lent);
% plot(f)
f=f+kpar(3)*sh_aif;
% figure;plot(f)

f=real(f);
if(~isfinite(f))
    f=zeros(1,lent);
end

