%
% 2-D cosine window
%
function wnd=cosw2d(sz);

w1=cosw1d(sz(1));
w2=cosw1d(sz(2));

wnd= w1(:)*w2(:)';

%
% 1-D cosine window
%
function wnd=cosw1d(nnm);

   nn=round(nnm/2);
   wnd=cosww(nn);
   if(nn*2==nnm)
     wnd=[flipud(wnd);wnd];
   else
     wnd=[flipud(wnd);wnd(2:end)];
   end


function w=cosww(n)
w=0.5+0.5*cos(pi*[0:n-1]'/(n-1));

