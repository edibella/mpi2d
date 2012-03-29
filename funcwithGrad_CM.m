% just here for eventual use of including gradients. Not yet for the mpi code

function [Err,Grad]=funcCM(x)


% x is a dummy for now
% matlab function to be called by constr

global nFactors
global avg
global ii
global N
global M
global measured_dixels
global truebld
global Fstart
global nDixels
global nFrames
global ttimes
global C



%disp('Iteration ')
ii=ii+1;

for i=1:nFactors
   C(1:nDixels,i) = x((i-1)*nDixels+1:i*nDixels);
end
kwo=x(nFactors*nDixels+1:nFactors*nDixels+nFactors-1)

% compute F from truebld and kparams:

F=zeros(nFrames,nFactors);
F(1:nFrames,1) = truebld;


sample_t=1;
for i=1:nFactors-1
  kernels(i,:) = (sample_t/60)*exp(-kwo(i)*sample_t/60*ttimes') ;
  tmp = conv(truebld,kernels(i,:));
  F(1:nFrames,i+1) = tmp(1:nFrames)';
end

dixels = C * F' ;


%disp('Sum of negatives in A and F is :')
negs_weight=20;
Fsum = -negs_weight*negs_weight*(sum(sum(C - abs(C))) + sum(sum(F - abs(F))));
Cnegs=C-abs(C);
Cnegs=Cnegs.*Cnegs;
%Fnegs=F-abs(F);
%Fnegs=Fnegs.*Fnegs;
Fsum = negs_weight*negs_weight*(sum(sum(Cnegs)));



%plot(dixels(5,:),'m')
%plot(measured_dixels(5,:),'b')
%pause
%plot(dixels(74,:),'mx')
%plot(measured_dixels(74,:),'bx')
%pause
Err = sum(sum( (dixels - measured_dixels) .* (dixels - measured_dixels) ));

% compute gradients of C and of kparams:
roughC=1;
e= dixels - measured_dixels; 
dfdC=zeros(nDixels,nFactors);
dfdk=zeros(1,nFactors-1);
for k=1:nFactors
   for i=1:nFrames
    	dfdC(:,k)=dfdC(:,k)+ 2* e(:,i).*F(i,k); 
   end
   dfdC(:,k)=dfdC(:,k)+ 2*roughC*roughC*C(:,k);  % non-neg. part
end

for  k=1:nFactors-1
   deriv_kernels(k,:) = -ttimes'.* kernels(k,:);
   tmpcurve=conv(truebld,deriv_kernels(k,:));
   tmpcurve=tmpcurve(1:nFrames)';
   for i=1:nFrames
      for j=1:nDixels
    	dfdk(k)=dfdk(k)+ 2* e(j,i).*C(j,k).*tmpcurve(i); 
      end
   end
end

% make grad a vector like x
for k=1:nFactors
   Grad(nDixels*(k-1)+1:nDixels*k) = dfdC(1:nDixels,k);
end
Grad(nDixels*nFactors+1:nDixels*nFactors+nFactors-1)=dfdk;

%e = (dixels - measured_dixels) .* (dixels - measured_dixels);
%disp('e is ')
%[I, J]=find(e>0);
%disp('I is '),I'
%disp('J is '),J'

%plot(F(1:nFrames,1),'b')
%plot(dixels(300,:),'k')
%plot(truedixels(300,:),'m')
%drawnow

Fsum=0;
Err = Err + Fsum

ii
%Err=x(1)*x(1)
%Err = Fsum

%if(ii%100 == 0)
%   disp('Iteration ')
%   ii
%   Err
%end


 
%disp('Max of max of abs of diff. of dixels from original is :')
%max(max(abs(dixels - truedixels)))
%disp('Sum of squared diff. of dixels from original is :')
%Fsum = sum(sum( (dixels - orig_dixels) .* (dixels - orig_dixels) ))
% add up all the negatives (they're doubled)...
% minimize the negative of that, so best possible is 0.

beta=10; % arbitrary weight introduced 9/25/00
%Fsum = (sum(sum(AA.*AA)) + beta*sum(sum(FF.*FF)) )
% sum of squares of F is 1
%Fsum = Fsum + 1 - sum(F(1,:).*F(1,:));
%I=find(F > 1.0 )
%Fsum = Fsum + sum(sum(F(I).*F(I)))

%I = find( F < 0 );
%F(I)=zeros(size(F(I)));
 

%disp('And quasi-rotation matrix B is :'), B
G=-10;
