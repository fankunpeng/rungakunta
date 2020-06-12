function [Cbp, chiCbp, Css, chiCss] = Chowtest(y, Tb, p, inc, bs)

%Cbp = Chow break point test statistic
%Css = Chow split sample test statistic

%y = (T+p)xK matrix of observations
%Tb = break date
%p = lag order
%inc = indicator for intercept
%bs = number of bootstrap replications


%Chow break point tests:

%points in time, def as 2-39
[Traw,K] = size(y);
rhs = lagmatrix(y,1:p);

% Y and Z for the whole sample:
Z12 = rhs(p+1:end,:)';
size(Z12)
T12 = length(Z12);
if inc == 1
    Z12 = [ones(1,T12); Z12];
end
Y12 = y(p+1:Traw,:)';

%  Y and Z for the  first subsample:  
Y1 = y(p+1:(Tb-1),:)'; % all obvs until Tb
T1 = length(Y1);
Y1
Z1 = rhs(p+1:(Tb-1),:)';
if inc == 1
    Z1 = [ones(1,T1); Z1];
end
Z1

 
%  Y and Z for the second subsample:  
Y2 = y(Tb+p:end,:)';  
T2 = length(Y2);
Z2 = rhs(Tb+p:Traw,:)';
if inc ==1
    Z2 = [ones(1,T2); Z2];
end
  
Btilde1 = Y1*Z1'*inv(Z1*Z1'); %LS Estimator
Btilde2 = Y2*Z2'*inv(Z2*Z2');
Btilde12 = Y12*Z12'*inv(Z12*Z12');


%set up the maximized log likelihood value of the unrestricted model

Utilde1 = Y1 - Btilde1*Z1;
Utilde2 = Y2 - Btilde2*Z2;
Uhat12 = (Y12 - Btilde12*Z12);

sigtilde1 = 1/T1*Utilde1*Utilde1';
sigtilde2 = 1/T2*Utilde2*Utilde2';

%loglM1 = -(T1/2)*log(abs(sigtilde1)) - (T2/2)*log(abs(sigtilde2));
loglM1 = -(T1/2)*log(det(sigtilde1)) - (T2/2)*log(det(sigtilde2)); % DETERMINANTE
 
%set up maximized log likelihood value of the restricted model (full sample)


%sigtilde12 = 1/(T1+T2) * ((Uhat12(:,p+1:T1)*Uhat12(:,p+1:T1)') + (Uhat12(:,T2:Traw-p)*Uhat12(:,T2:Traw-p)'));
sigtilde12 = 1./(T1+T2)*( Uhat12(:,1:T1)*Uhat12(:,1:T1)' + Uhat12(:,(T12-T2+1):T12)*Uhat12(:,(T12-T2+1):T12)' ) ;


loglM2 = -(T1+T2)/2*log(det(sigtilde12)); %DETERMINANTE

%Chow break point test statistic

Cbp = 2*(loglM1 - loglM2);

dofCbp = K^2*p+K+(K*(K+1))/2

chiCbp = 1 - chi2cdf(Cbp,dofCbp); % P VALUE 
%chi2inv(0.95, dofCbp);

%Chow sample split test

%loglM3 = -(T1+T2)/2 * log(det((1/(T1+T2))*(T1*sigtilde1 + T2*sigtilde2)));

Css = (T1 + T2) * (log(det(sigtilde12)) - log(det((T1+T2)^-1*(T1*sigtilde1 + T2*sigtilde2))));

dofCss = K^2*p+K;

chiCss =  1 - chi2cdf(Css,dofCss); % P VALUE
%chi2inv(0.95,dofCss);






end

