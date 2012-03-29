function [COEF,Serr,XTXI,Rsq,Fval,TTAB,FITS,RES] = lregress(Y,X)

% LREGRESS  Performs multiple linear regression analysis of Y on X.
%
%           Usage:
%
%           [ COEF, Serr, XTXI, Rsq, Fval, TTAB, FITS, RES ] = lregress(Y,X)
%
%           Returns:
%
%           COEF - Coefficients
%           Serr - Standard error of estimate
%           XTXI - inverse of X'*X
%           Rsq  - R-squared
%           Fval - F-value
%           TTAB - Coefficients, standard deviations, T-values
%           FITS - Fitted values
%           RES  - Residuals

% M J Chlond - Nov 93
% m.chlond@uclan.ac.uk


%  identify dimensions

    DIM = size(X) ;
    n = DIM(1,1)  ;
    k = DIM(1,2)  ;

%  solve normal equations

    X = [ ones(n,1) X ]   ;
    XTXI = inv(X'*X);
    COEF = XTXI*X'*Y ;

%  fitted values

    FITS = X*COEF ;

%  R-squared

    my = sum(Y)/n             ;
    TL = FITS' - ones(1,n)*my ;
    tl = TL*TL'               ;
    BL = Y - ones(n,1)*my     ;
    bl = BL'*BL               ;
    Rsq = tl/bl               ;

%  residuals and standard error

    RES = Y - FITS               ;
    Serr = (RES'*RES/(n-k-1))^.5 ;

%  coeffs, s.d. of slopes and t-values

    C = Serr^2*inv(X'*X)         ;
    C = sqrt(diag(C,0))       ;
    TTAB = [ COEF,C,COEF./C ] ;

%  F-value

    Fval = (tl/k)*(n-k-1)/(bl-tl) ;

%  end of procedure
