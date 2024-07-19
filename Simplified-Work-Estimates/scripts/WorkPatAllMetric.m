numPatients = 8;
sname = {'BiV1','BiV2','BiV3','BiV4','BiV5','BiV6','BiV7','BiV8'};

PatientData = load('PatientData.txt');
% PatientData = PatientData(:,[1:5,7:8]);
fileID = fopen('PatientDataLabels.txt');
PatientDataLabels = textscan(fileID,'%q');
PatientDataLabels = PatientDataLabels';
fclose(fileID);
currentsname = sname;

followUpDays = PatientData(15,:);
figure;
for plotNum = 4%1:14
%  	subplot(4,4,plotNum)
	patMetric = 7;%
% 	currentMetric = PatientData(patMetric,:)./followUpDays;
	currentMetric = 100.*PatientData(patMetric,:);
% 	currentFraction = covWork(1,:);
 	currentFraction = cell2mat(lvNegFraction);
 	currentFraction = currentFraction(1:8);
%  	currentFraction = sepNegFractionImp;
% 	currentFraction = (currentFraction(1:8)-currentFraction(9:16))./(currentFraction(1:8));
% 	currentFraction = covWork(1,1:8);
% 	currentFraction = (sumNegWorkImpFrac(1:7));
% 	xLabel = 'Septum Negative Work Volume Fraction';
% 	xLabel = 'COV Work';
	xLabel = workMetricLabels{2}; %15
	yLabel = PatientDataLabels{1}(patMetric);
% 	currentMetric(5) = [];
% 	currentFraction(5) = [];
% 	currentsname(5) = [];

	PlotCorrellationPaper(currentFraction,currentMetric,xLabel,yLabel,currentsname);
end