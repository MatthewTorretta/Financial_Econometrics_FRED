clc; clear;


%read in data from an excel file
[data,names]=xlsread('data.xls');

time=data(:,1);
Export=data(:,2);
Import=data(:,3);
ts1=[Export,Import];
plot(time, ts1)
legend('Export','Import')


Y=Export(1:130,:);

T=100;

ToEstMdlm1 = arima('ARLags',1,'MALags',1);
[EstMdlm1, logL1] = estimate(ToEstMdlm1,Y(1:T));



h = 30; %%%horizons
[YFm1,YMSEm1] = forecast(EstMdlm1,h,'Y0',Y(1:T));

%%%% calculate forecast error
EFm1 = Y(T+1:T+h)-YFm1;

%%%% calculate Root Mean Square Forecast Error

EF2m1=EFm1.^2;

RMSFEm1 = [sum(EF2m1)]/h;

figure
h1 = plot(Y,'Color',[.7,.7,.7]);
hold on
h2 = plot(101:101+h-1,YFm1,'b','LineWidth',2);
h3 = plot(101:101+h-1,YFm1 + 1.96*sqrt(YMSEm1),'r:',...
		'LineWidth',2);
plot(101:101+h-1,YFm1 - 1.96*sqrt(YMSEm1),'r:','LineWidth',2);
legend([h1 h2 h3],'Observed','Forecast',...
		'95% Confidence Interval','Location','NorthWest');
title(['30-Period Forecasts and Approximate 95% '...
			'Confidence Intervals'])
hold off




ToEstMdlm2 = arima('ARLags',1:2,'MALags',1:2);
EstMdlm2 = estimate(ToEstMdlm2,Y(1:T));

h = 30; %%%horizons
[YFm2,YMSEm2] = forecast(EstMdlm2,h,'Y0',Y(1:T));

%%%% calculate forecast error
EFm2 = Y(T+1:T+h)-YFm2;

%%%% calculate Root Mean Square Forecast Error

EF2m2=EFm2.^2;

RMSFEm2 = [sum(EF2m2)]/h;

figure
h1 = plot(Y,'Color',[.7,.7,.7]);
hold on
h2 = plot(101:101+h-1,YFm2,'b','LineWidth',2);
h3 = plot(101:101+h-1,YFm2 + 1.96*sqrt(YMSEm2),'r:',...
		'LineWidth',2);
plot(101:101+h-1,YFm2 - 1.96*sqrt(YMSEm2),'r:','LineWidth',2);
legend([h1 h2 h3],'Observed','Forecast',...
		'95% Confidence Interval','Location','NorthWest');
title(['30-Period Forecasts and Approximate 95% '...
			'Confidence Intervals'])
hold off


%[stat, prob]  = DMtest(EFm1, EFm2, h)

 [DM, prob] = dmtest1(EFm1, EFm2, h)