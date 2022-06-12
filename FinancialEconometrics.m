clear;
clc;

cd

addpath(genpath(strcat(cd,'\data')));    
addpath(genpath(strcat(cd,'\functions'))); 
addpath(genpath(strcat(cd,'\jplv7')));

%% Import an excel file

[data,text]= xlsread ('fred_qd.xlsx');
names=text(1,:);

% FRED_QD dataset contains macroeconomic sampled at a quarterly frequency
% Data are from 1959-Q1 to 2019-Q4

%% Select the sample 
DATAFINAL=data;

%% Select variables of interest (i.e. 
% 1. PCECTPI:  Personal Consumption Expenditures: Chain-type Price Index (Index 2012=100)  
% 2. DFXARG3Q086SBEA: Personal consumption expenditures: Nondurable goods: Food and beverages purchased for off-premises consumption (chain-type price index) 
% 3. DHLCRG3Q086SBEA: Personal consumption expenditures: Services: Health care (chain-type price index)

labx=["PCECTPI" "DFXARG3Q086SBEA" "DHLCRG3Q086SBEA"];

idvars=find(contains(names,labx));

                                   
DATA=DATAFINAL(:,idvars);

%% a) OLS estimation.
% The model can be written as follows
% PCE = b0 + b1*FOOD&BEV + b2*HEALTHCARE + epsilon
% or in matrix notation:
%   Y   =   X    *   beta    +  EPSILON
% (T*1)   T*(N+1)  (N+1)*1       (T*1)

Y=DATA(:,1);       % dependent variables
X=DATA(:,2:end);   % independent variables

alpha=0.05; % significance level of the test

[T,N]=size(X);

X=[ones(T,1) X];  
ncoef=size(X,2);  
dof=T-ncoef;   

beta= inv(X'*X)*X'*Y % beta OLS

epsilon=Y-X*beta; 
residsummary=[min(epsilon) quantile(epsilon,[.25,.5,.75]) max(epsilon)];
%% b) Confidence Interval.

RSS=(epsilon'*epsilon);    
sigmaeps=(1./dof).*RSS;   
rse=sqrt(sigmaeps);

varbeta=sigmaeps.*inv(X'*X); 
stdbeta=sqrt(diag(varbeta)); 

tvals=beta./stdbeta; %t-statistics
tcrit=-1.*tinv(alpha./2,dof);

confint=[beta-tcrit.*stdbeta beta+tcrit.*stdbeta]  % Confidence intervals

%% c) Statistical Significance.

pvals=2.*(1-tcdf(abs(tvals),dof)) % compute p-values

%% d) Rsquared and Adjusted Rsquared.

YDEMEANED=Y-mean(Y); 
TSS=YDEMEANED'*YDEMEANED; 

R2=1-(RSS./TSS) % R-squared

RSSBAR=RSS./dof;
TSSBAR=TSS./(T-1);

R2BAR=1-(RSSBAR./TSSBAR) % Adjusted R-squared


%% e) F-statistic.

% Restricted model 
XTILDE=ones(T,1);
betatilde=inv(XTILDE'*XTILDE)*XTILDE'*Y; 
epstilde=Y-XTILDE*betatilde;
RSSRESTR=epstilde'*epstilde;  
RSSUNRESTR=RSS;

% F-statistic
Fstat=((RSSRESTR-RSSUNRESTR)./(ncoef-1))./(RSSUNRESTR./(dof))
Fpval=1-fcdf(Fstat,ncoef-1,dof)

%% f) Plot of the fitted model.

figure(2);
plot(epsilon)
hold on
axis tight

% Fitted vs. residuals
figure(3);
yfitted=X*beta;
scatter(yfitted,epsilon,'filled')
hold on;
plot(sort(yfitted),zeros(size(epsilon,1),1),'b--');
axis tight
title('Residuals vs. Fitted');
xlabel('Fitted'); ylabel('Residuals')

%% g) Diagnostic tests for residuals (normality).

%Normality
figure(4)
subplot(1,2,1)
qqplot(epsilon)
title('Normal Q-Q');
ylabel('Residuals');
subplot(1,2,2)
histogram(epsilon);
title('Histogram of Residuals');
xlabel('Residuals'); ylabel('Frequency');

% Jarque-Bera test
alpha=0.05; 
[h,p,jbstat,critval]=jbtest(epsilon,alpha)

%% g) Diagnostic tests for residuals (Heteroskedacity).

RSS=(epsilon'*epsilon);   
                          
sigmaeps=(1./T).*RSS;     
                          

epsnew=(epsilon.^2)./sigmaeps - 1;  
Z=X
betatilde=inv(Z'*Z)*Z'*epsnew; 
epsfitted=Z*betatilde; 
BP=sum(epsfitted.^2)./2 

bpcrit=chi2inv(0.95,N) 
pvalue=1-chis_prb(BP,N) 

bpagan(Y,X)

%% g) Diagnostic tests for residuals (Homoskedacity).

ZTILDE=X;
for i=1:N
    for j=1:N
       if i==j
          ZTILDE=[ZTILDE X(:,i+1).^2];
       elseif i<j
          ZTILDE=[ZTILDE X(:,i+1).*X(:,j+1)];
       end 
    end
end

betanew=inv(ZTILDE'*ZTILDE)*ZTILDE'*(epsilon.^2); 
uresid=(epsilon.^2)-ZTILDE*betanew;

YDEMEANED=(epsilon.^2)-mean((epsilon.^2));
TSS=YDEMEANED'*YDEMEANED; 
RSSNEW=(uresid'*uresid);

R2=1-(RSSNEW./TSS); 
W=T.*R2

Wcrit=chi2inv(0.95,size(ZTILDE,2)) 

%% g) Diagnostic tests for residuals (Serial Correlation).

DW=sum(diff(epsilon,1).^2)./sum(epsilon.^2)

%% h) Multicollinearity.

XAUGMENTED=[X yfitted.^2 yfitted.^3]; 
[TSTAR,NSTAR]=size(XAUGMENTED);

BETAAUGMENTED=inv(XAUGMENTED'*XAUGMENTED)*XAUGMENTED'*Y; 
epsaugmented=Y-XAUGMENTED*BETAAUGMENTED; 

RSSUNRESTR=epsaugmented'*epsaugmented;  

RSSRESTR=RSS; 

jrestriction= size(XAUGMENTED,2)- size(X,2);  

dofstar=TSTAR-NSTAR;

Fstat=((RSSRESTR-RSSUNRESTR)./(jrestriction))./(RSSUNRESTR./(dofstar))
Fpval=1-fcdf(Fstat,jrestriction,dofstar)

