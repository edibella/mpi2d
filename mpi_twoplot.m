function [coeff] = mpi_twoplot(Y,X,offset)

% twoplot    Plot regression of Y on X and X on Y.
%
%            Usage:
%
%            [coeff] = twoplot(Y,X)
%
%            Returns:
%
%            coeff  - column 1 gives coefficients of regression Y on X
%                     column 2 gives rearranged coefficients of X on Y

% M J Chlond - May 96
% m.chlond@uclan.ac.uk

% modified a little to print out t-values of slope as well
%coeffyx = lregress(Y,X);
[COEF,Serr,XTXI,Rsq,Fval,TTAB,FITS,RES]= lregress(Y,X);

coeffyx = COEF;
%disp('t value is (bottom row is for slope): ')
%TTAB(:,3)
disp('coeffs., s.d. of slopes, t values (bottom row is for slope): ')
TTAB
disp('Serr is '), Serr

coeffxy = lregress(X,Y);
coeffxy = [-coeffxy(1)/coeffxy(2);1/coeffxy(2)]; 
coeff = [coeffyx coeffxy];

minx = min(X);
maxx = max(X)+4;
t1 = linspace(minx,maxx,16);
t1 = linspace(0,maxx,16);
XI = [ones(1,16);t1];
YI = coeff'*XI;

%plot(X+offset-1,Y,'+',XI(2,:)+offset-1,YI(1,:),'linewidth',2.5)
plot(XI(2,:)+offset-1,YI(1,:),'k','linewidth',3)
% hmmm, line seems a little crooked at top...
% evd changed to just plot out regression of Y on X
%plot(X,Y,'+',XI(2,:),YI(1,:),XI(2,:),YI(2,:),':')

%%%% to automatically display the linear regression eqns and R^2 values on the plot 
xx=3;yy=1.2; %location of eqn displayed on plot
labelString1=sprintf('y=%1.2f*x+%1.2f',coeffyx(2),coeffyx(1));
text(xx,yy,labelString1,'fontsize',16,'fontweight','bold')
xx=3;yy=0.6;
labelString2=sprintf('r=%1.2f',sqrt(Rsq));
text(xx,yy,labelString2,'fontsize',16,'fontweight','bold')

disp('Rsq is ')
Rsq

disp('R is ')
sqrt(Rsq)

disp('Fval is ')
Fval

return
