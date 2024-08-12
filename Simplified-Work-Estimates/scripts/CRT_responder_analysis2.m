% This script evaluates the ability of our simple MW estimates to predict
% CRT response. This is insipred by Adarsh's analysis
% 1. work
% 2. coefficient of variation (COV) of work
% 3. LV fraction performing negative work
% 4. LV septal fraction perofrming negative work?

% lines 314 - 347 plot the figures included in the manuscript

%to add: fishers transformation to evaluate how different the relatinships
%are to the ground truth
clear; clc


% NOTE: make sure your home path is the Simplified Work Estimates directory
homepath = '/Users/amandacraine/Documents/ContijochLab/repos/PSM-Work-Mapping/Simplified-Work-Estimates';
addpath([homepath,'/data/'])
addpath([homepath,'/scripts'])
addpath([homepath,'/LV Geometric Models'])
addpath([homepath,'/figures'])
figpath = [homepath,'/figures/'];
addpath([homepath,'/scripts/CRT Analysis Scripts/'])
cd(homepath)

%need the work values for all points
load("WorkPSMLBBBCRT.mat")

%PSA work 
load('MWCT_tris.mat')
load('Laplace_work_all_patches.mat')
load('GenericLVP_work_all_patches.mat')
load('PatchAreas.mat')
%load('segMWCT_effrad_allpats.mat')
%load('all_seg_work_all_pats.mat')
figtitle = [{'workPSM'},{'MWCT'},{'MWCTLP'},{'MWCTLPTV'},{'MWCTgenLVP'},{'MWCTgenLVP_ESP'}];

%%
for pat = 1:8
    workPSM{pat} = work{pat}(gauss{1}).*detj{pat}(gauss{1});

    %total work (in J or Pa)
    tot_workPSM(pat) = sum(workPSM{pat})./1000;
    tot_workMWCT(pat) = sum(MWCT_tris{pat})./1000;
    tot_workMWCTLP(pat) = sum(MWCTLP(:,pat))./1000;
    tot_workMWCTLPTV(pat) = sum(MWCTLPTV(:,pat))./1000;
    tot_workMWCTgenericLVP(pat) = sum(MWCTgenericLVP(:,pat))./1000;
    tot_workMWCTgenericLVP_patES(pat) = sum(MWCTgenericLVP_patES(:,pat))./1000;

    %coefficient of variation regional work (COVW)
    cov_workPSM(pat) = stdWork(pat,1)./meanWork(pat,1);%std(workPSM{pat})./mean(workPSM{pat});
    cov_workMWCT(pat) = std(MWCT_tris{pat})./mean(MWCT_tris{pat});
    cov_workMWCTLP(pat) = std(MWCTLP(:,pat))./mean(MWCTLP(:,pat));
    cov_workMWCTLPTV(pat) = std(MWCTLPTV(:,pat))./mean(MWCTLPTV(:,pat));
    cov_workMWCT_genericLVP(pat) = std(MWCTgenericLVP(:,pat))./mean(MWCTgenericLVP(:,pat));
    cov_workMWCT_genericLVP_patES(pat) = std(MWCTgenericLVP_patES(:,pat))./mean(MWCTgenericLVP_patES(:,pat));
    
    meanMWCT{pat} = sum(MWCT_tris{pat}'.*PatchAreas{pat}(:,1))/sum(PatchAreas{pat}(:,1));
    solnMWCT = MWCT_tris{pat}' - meanMWCT{pat};
    stdMWCT{pat} = sqrt(sum((solnMWCT.^2).*PatchAreas{pat}(:,1))/sum(PatchAreas{pat}(:,1)));

    meanMWCTLP{pat} = sum(MWCTLP(:,pat).*PatchAreas{pat}(:,1))/sum(PatchAreas{pat}(:,1));
    solnMWCTLP = MWCTLP(:,pat) - meanMWCTLP{pat};
    stdMWCTLP{pat} = sqrt(sum((solnMWCTLP.^2).*PatchAreas{pat}(:,1))/sum(PatchAreas{pat}(:,1)));

    meanMWCTLPTV{pat} = sum(MWCTLPTV(:,pat).*PatchAreas{pat}(:,1))/sum(PatchAreas{pat}(:,1));
    solnMWCTLPTV = MWCTLPTV(:,pat) - meanMWCTLPTV{pat};
    stdMWCTLPTV{pat} = sqrt(sum((solnMWCTLPTV.^2).*PatchAreas{pat}(:,1))/sum(PatchAreas{pat}(:,1)));

    meanMWCTgenericLVP{pat} = sum(MWCTgenericLVP(:,pat).*PatchAreas{pat}(:,1))/sum(PatchAreas{pat}(:,1));
    solnMWCTgenericLVP = MWCTgenericLVP(:,pat) - meanMWCTgenericLVP{pat};
    stdMWCTgenericLVP{pat} = sqrt(sum((solnMWCTgenericLVP.^2).*PatchAreas{pat}(:,1))/sum(PatchAreas{pat}(:,1)));

    meanMWCTgenericLVP_patES{pat} = sum(MWCTgenericLVP_patES(:,pat).*PatchAreas{pat}(:,1))/sum(PatchAreas{pat}(:,1));
    solnMWCTgenericLVP_patES = MWCTgenericLVP_patES(:,pat) - meanMWCTgenericLVP_patES{pat};
    stdMWCTgenericLVP_patES{pat} = sqrt(sum((solnMWCTgenericLVP_patES.^2).*PatchAreas{pat}(:,1))/sum(PatchAreas{pat}(:,1)));

    covMWCT(pat) = stdMWCT{pat}./meanMWCT{pat};
    covMWCTLP(pat) = stdMWCTLP{pat}./meanMWCTLP{pat};
    covMWCTLPTV(pat) = stdMWCTLPTV{pat}./meanMWCTLPTV{pat};
    covMWCTgenericLVP(pat) = stdMWCTgenericLVP{pat}./meanMWCTgenericLVP{pat};
    covMWCTgenericLVP_patES(pat) = stdMWCTgenericLVP_patES{pat}./meanMWCTgenericLVP_patES{pat};

    %LV negative work
    lvVolume			= detj{pat}(gauss{1});
    lvWork				= work{pat}(gauss{1});
    lvNegLocation		= find(lvWork < 0.0);
    lvNegVolume			= sum(lvVolume(lvNegLocation));
	lvNegFractionPSM(pat)	= lvNegVolume/sum(lvVolume);

    %lvMWCTNegLocation   = find(MWCT_tris{pat} < 0);
    %lvMWCTNegVolume     = sum(PatchAreas{pat}(lvMWCTNegLocation,1));
    %lvNegFractionMWCT(pat) = lvMWCTNegVolume./sum(PatchAreas{pat}(:,1));
    lvNegFractionMWCT(pat) = getNegWorkFractions(MWCT_tris{pat},PatchAreas{pat}(:,1));
    lvNegFractionMWCTLP(pat) = getNegWorkFractions(MWCTLP(:,pat),PatchAreas{pat}(:,1));
    lvNegFractionMWCTLPTV(pat) = getNegWorkFractions(MWCTLPTV(:,pat),PatchAreas{pat}(:,1));
    lvNegFractionMWCT_genericLVP(pat) = getNegWorkFractions(MWCTgenericLVP(:,pat),PatchAreas{pat}(:,1));
    lvNegFractionMWCT_genericLVP_patES(pat) = getNegWorkFractions(MWCTgenericLVP_patES(:,pat),PatchAreas{pat}(:,1));

    % LVnegworkPSM(pat) = numel(find(workPSM{pat} < 0))./numel(workPSM{pat});
    % LVnegworkMWCT(pat) = numel(find(MWCT_tris{pat} < 0))./numel(MWCT_tris{pat});
    % LVnegworkMWCTLP(pat) = numel(find(MWCTLP(:,pat) < 0))./numel(MWCTLP(:,pat));
    % LVnegworkMWCTLPTV(pat) = numel(find(MWCTLPTV(:,pat) < 0))./numel(MWCTLPTV(:,pat));
    % LVnegworkMWCT_genericLVP(pat) = numel(find(MWCTgenericLVP(:,pat) < 0))./numel(MWCTgenericLVP(:,pat));
    % LVnegworkMWCT_genericLVP_patES(pat) = numel(find(MWCTgenericLVP_patES(:,pat) < 0))./numel(MWCTgenericLVP_patES(:,pat));

    %septum negative work
    sepVolume			= detj{pat}(gauss{3});
	sepWork				= work{pat}(gauss{3});
	sepNegLocation		= find(sepWork < 0.0);
	sepNegVolume		= sum(sepVolume(sepNegLocation));
	sepNegFractionPSM(pat)	= sepNegVolume/sum(sepVolume);

    [pointsSeg,TriSegElems] = getSegmentLabels(elementIndices,pat,1);
    sep_indx = find(TriSegElems == 2 |TriSegElems == 3 |TriSegElems == 8 |TriSegElems == 9 |TriSegElems == 14);
    
    sepNegFractionMWCT(pat) = getNegWorkFractions(MWCT_tris{pat}(sep_indx),PatchAreas{pat}(sep_indx,1));
    sepNegFractionMWCTLP(pat) = getNegWorkFractions(MWCTLP(sep_indx,pat),PatchAreas{pat}(sep_indx,1));
    sepNegFractionMWCTLPTV(pat) = getNegWorkFractions(MWCTLPTV(sep_indx,pat),PatchAreas{pat}(sep_indx,1));
    sepNegFractionMWCT_genericLVP(pat) = getNegWorkFractions(MWCTgenericLVP(sep_indx,pat),PatchAreas{pat}(sep_indx,1));
    sepNegFractionMWCT_genericLVP_patES(pat) = getNegWorkFractions(MWCTgenericLVP_patES(sep_indx,pat),PatchAreas{pat}(sep_indx,1));
   
    % segsegworkPSM(pat) = numel(find(work{pat}(gauss{3}).*detj{pat}(gauss{3})<0))./numel(work{pat}(gauss{3}).*detj{pat}(gauss{3}));
    % sepnegworkMWCT(pat) = numel(find(MWCT_tris{pat}(sep_indx) < 0))./numel(MWCT_tris{pat}(sep_indx));
    % sepnegworkMWCTLP(pat) = numel(find(MWCTLP(sep_indx,pat) < 0))./numel(MWCTLP(sep_indx,pat));
    % sepnegworkMWCTLPTV(pat) = numel(find(MWCTLPTV(sep_indx,pat) < 0))./numel(MWCTLPTV(sep_indx,pat));
    % sepnegworkMWCT_genericLVP(pat) = numel(find(MWCTgenericLVP(sep_indx,pat) < 0))./numel(MWCTgenericLVP(sep_indx,pat));
    % sepnegworkMWCT_genericLVP_patES(pat) = numel(find(MWCTgenericLVP_patES(sep_indx,pat) < 0))./numel(MWCTgenericLVP_patES(sep_indx,pat));
end
%% total work
totwork = [tot_workPSM;tot_workMWCT;tot_workMWCTLP;tot_workMWCTLPTV;...
    tot_workMWCTgenericLVP;tot_workMWCTgenericLVP_patES];
ymax = max(totwork,[],"all") + 0.15.*range(totwork,"all");
ymin = min(totwork,[],"all") - 0.15.*range(totwork,"all");
titles = [{'SSA_P_S_M'},{'P_L_H_CSA'},{'WS_E_DSA'},{'WS_T_VSA'},...
    {'P_g_e_n_e_r_i_cSA'},{'P_g_e_n_e_r_i_c_,_s_c_a_l_e_dSA'}];
figtitle = [{'workPSM'},{'MWCT'},{'MWCTLP'},{'MWCTLPTV'},{'MWCTgenLVP'},{'MWCTgenLVP_ESP'}];
ytitle = 'LV Total Myocardial Work';
for i = 1:length(titles)
    % if i == 1
    %     ytitle = 'Myocardial Work (J)';
    % else
    %     ytitle = 'Myocardial Work (Pa)';
    % end
    %fig = plotResponderBoxplot(totwork(i,:),'LV',titles{i},ytitle,ymin,ymax);
    %saveas(fig,['CRT Analysis Plots/total_work_',figtitle{i},'.png'])
end
fig = plotResponderBoxplotAllEsimates(totwork(1,:),totwork(2:end,:),ymin,ymax,ytitle);

%%
totwork = [tot_workPSM;tot_workMWCT;tot_workMWCTLP;tot_workMWCTLPTV;...
    tot_workMWCTgenericLVP;tot_workMWCTgenericLVP_patES];
ymax = max(totwork,[],"all") + 0.15.*range(totwork,"all");
ymin = min(totwork,[],"all") - 0.15.*range(totwork,"all");
xmin = min(tot_workPSM) - 0.15.*range(tot_workPSM);
xmax = max(tot_workPSM) + 0.15.*range(tot_workPSM);
titles = [{'SSA_P_S_M'},{'P_L_H_CSA'},{'WS_E_DSA'},{'WS_T_VSA'},...
    {'P_g_e_n_e_r_i_cSA'},{'P_g_e_n_e_r_i_c_,_s_c_a_l_e_dSA'}];
figtitle = [{'workPSM'},{'MWCT'},{'MWCTLP'},{'MWCTLPTV'},{'MWCTgenLVP'},{'MWCTgenLVP_ESP'}];
ytitle = [{''},{'Total P_L_H_CSA (Pa)'},{'Total WS_E_DSA (Pa)'},{'Total WS_T_VSA (Pa)'},...
    {'Total P_g_e_n_e_r_i_cSA (Pa)'},{'Total P_g_e_n_e_r_i_c_,_s_c_a_l_e_dSA (Pa)'}];
xtitle = 'Total PSM Myocardial Work (J)';
for i = 2:size(totwork,1)
    fig = plotCorrelationPatientData(totwork(1,:),totwork(i,:),titles(i),xtitle,ytitle(i),ymin,ymax,xmin,xmax);
    %saveas(fig,['CRT Analysis Plots/totwork_corr_',figtitle{i},'.png'])
end

% figure; hold all
% plot(totwork(1,:),totwork(2,:),'.','MarkerSize',20)
% plot([0 2],[0 2],'k')
% for i = 1:6
%     r(i) = corr(totwork(1,:)',totwork(i,:)');
% end
%% COV
%fig1 = plotResponderBoxplotCOV(cov_workPSM,cov_workCT,'LV',0,1);
location = 'LV';
ytitle = 'COV';
figtitle = [{'workPSM'},{'MWCT'},{'MWCTLP'},{'MWCTLPTV'},{'MWCTgenLVP'},{'MWCTgenLVP_ESP'}];
cov = [cov_workPSM;cov_workMWCT;cov_workMWCTLP;cov_workMWCTLPTV;...
    cov_workMWCT_genericLVP;cov_workMWCT_genericLVP_patES];
titles = [{'SSA_P_S_M'},{'P_L_H_CSA'},{'WS_E_DSA'},{'WS_T_VSA'},...
    {'P_g_e_n_e_r_i_cSA'},{'P_g_e_n_e_r_i_c_,_s_c_a_l_e_dSA'}];
ymax = max(cov,[],"all") + 0.15.*range(cov,"all");
ymin = min(cov,[],"all") - 0.15.*range(cov,"all");
xmin = ymin;%min(cov(1,:)) - 0.15*range(cov(1,:));
xmax = ymax;%max(cov(1,:)) + 0.15*range(cov(1,:));
% for i = 1%1:length(titles)
%     fig = plotResponderBoxplot(cov(i,:),'LV',titles{i},ytitle,ymin,ymax);
%     %saveas(fig,['CRT Analysis Plots/cov_',figtitle{i},'.png'])
% end
% 
ytitle = [{''},{'COV P_L_H_CSA'},{'COV WS_E_DSA'},{'COV WS_T_VSA'},...
    {'COV P_g_e_n_e_r_i_cSA'},{'COV P_g_e_n_e_r_i_c_,_s_c_a_l_e_dSA'}];
xtitle = 'COV SSA_P_S_M';
for i = 2:size(cov,1)
    fig = plotCorrelationPatientData(cov(1,:),cov(i,:),titles(i),xtitle,ytitle(i),ymin,ymax,xmin,xmax);
    %saveas(fig,['CRT Analysis Plots/cov_corr_',figtitle{i},'.png'])
end

%% LV negative work
%fig2 = plotResponderBoxplot(lvNegFraction,LVnegworkCT,'LV',1,0);

LVnegwork = [lvNegFractionPSM;lvNegFractionMWCT;LVnegworkMWCT;LVnegworkMWCTLP;LVnegworkMWCTLPTV;...
    LVnegworkMWCT_genericLVP;LVnegworkMWCT_genericLVP_patES];
titles = [{'SSA_P_S_M'},{'V_f_LVNW MWCT'},{'P_L_H_CSA'},{'WS_E_DSA'},{'WS_T_VSA'},...
    {'P_g_e_n_e_r_i_cSA'},{'P_g_e_n_e_r_i_c_,_s_c_a_l_e_dSA'}];
ymax = max(LVnegwork,[],"all") + 0.15.*range(LVnegwork,"all");
ymin = min(LVnegwork,[],"all") - 0.15.*range(LVnegwork,"all");
xmin = ymin;
xmax = ymax;
% for i = 1:length(titles)
%     if i == 1
%         ytitle = 'LV Negative Work Volume Fraction';
%     else
%          ytitle = 'LV Negative Work Fraction';
%     end
%     fig = plotResponderBoxplot(LVnegwork(i,:),'LV',titles{i},ytitle,ymin,ymax);
%     saveas(fig,['CRT Analysis Plots/LVnegwork_',figtitle{i},'.png'])
% end

ytitle = [{''},{'LV Negative P_L_H_CSA Fraction'},...
    {'LV Negative WS_E_DSA Fraction'},{'LV Negative WS_T_VSA Fraction'},...
    {'LV Negative P_g_e_n_e_r_i_cSA Fraction'},{'LV Negative P_g_e_n_e_r_i_c_,_s_c_a_l_e_dSA Fraction'}];
xtitle = 'V_fLVNW LBBB';
for i = 2:size(LVnegwork,1)
    fig = plotCorrelationPatientData(LVnegwork(1,:),LVnegwork(i,:),titles(i),xtitle,ytitle(i),ymin,ymax,xmin,xmax);
   % saveas(fig,['CRT Analysis Plots/lvnegwork_corr_',figtitle{i},'.png'])
end

%% Septum negative work
%fig3 = plotResponderBoxplot(sepNegFraction,sepnegworkCT,'Septum',1,0);
sepnegwork = [sepNegFractionPSM;segsegworkPSM;sepnegworkMWCT;sepnegworkMWCTLP;sepnegworkMWCTLPTV;...
    sepnegworkMWCT_genericLVP;sepnegworkMWCT_genericLVP_patES];
titles = [{'SSA_P_S_M'},{'SSA_P_S_M'},{'P_L_H_CSA'},{'WS_E_DSA'},{'WS_T_VSA'},...
    {'P_g_e_nSA'},{'P_g_e_n_,_s_c_a_l_e_dSA'}];
ymax = max(sepnegwork,[],"all") + 0.15.*range(sepnegwork,"all");
ymin = min(sepnegwork,[],"all") - 0.15.*range(sepnegwork,"all");
xmin = ymin; xmax = ymax;
% for i = 1:length(titles)
%     if i == 1
%         ytitle = 'Septal Negative Work Volume Fraction';
%     else
%          ytitle = 'Fraction Septal Negative Work Fraction';
%     end
%     fig = plotResponderBoxplot(sepnegwork(i,:),'LV',titles{i},ytitle,ymin,ymax);
%     saveas(fig,['CRT Analysis Plots/sepnegwork_',figtitle{i},'.png'])
% end

ytitle = [{''},{'Septal Negative P_L_H_CSA Fraction'},...
    {'Septal Negative WS_E_DSA Fraction'},{'Septal Negative WS_T_VSA Fraction'},...
    {'Septal Negative P_g_e_nSA Fraction'},{'Septal Negative P_g_e_n_,_s_c_a_l_e_dSA Fraction'}];
xtitle = 'V_fSTNW LBBB';
for i = 2:size(sepnegwork,1)
    fig = plotCorrelationPatientData(sepnegwork(1,:),sepnegwork(i,:),titles(i),xtitle,ytitle(i),ymin,ymax,xmin,xmax);
    %saveas(fig,['CRT Analysis Plots/sepnegwork_corr_',figtitle{i},'.png'])
end

% for i = 1:6
%     r(i) = corr(sepnegwork(1,:)',sepnegwork(i,:)');
% end

%%
cd('data/')
load("WorkPSMLBBBCRT.mat")
cd ('../scripts')
WorkMetrics
WorkPatAllMetric

%% correlation between negative work and change in ESV
meas_idx = 1; %1 for LV, 2 for septum

% LVnegwork = [lvNegFractionPSM;LVnegworkMWCT;LVnegworkMWCTLP;LVnegworkMWCTLPTV;...
%     LVnegworkMWCT_genericLVP;LVnegworkMWCT_genericLVP_patES];
LVnegwork = [lvNegFractionPSM;lvNegFractionMWCT;lvNegFractionMWCTLP;lvNegFractionMWCTLPTV;...
    lvNegFractionMWCT_genericLVP;lvNegFractionMWCT_genericLVP_patES];
% sepnegwork = [sepNegFractionPSM;sepnegworkMWCT;sepnegworkMWCTLP;sepnegworkMWCTLPTV;...
%     sepnegworkMWCT_genericLVP;sepnegworkMWCT_genericLVP_patES];
sepnegwork = [sepNegFractionPSM;sepNegFractionMWCT;sepNegFractionMWCTLP;sepNegFractionMWCTLPTV;...
    sepNegFractionMWCT_genericLVP;sepNegFractionMWCT_genericLVP_patES];
patMetric = 7;

if meas_idx == 1
    location = 'LV';
    workmetric = LVnegwork;
elseif meas_idx == 2
    location = 'Septal'; %'LV'
    workmetric = sepnegwork;%LVnegwork;
end

deltaESV = PatientData(patMetric,:);
titles = [{'SSA_P_S_M'},{'P_L_H_CSA'},{'WS_E_DSA'},{'WS_T_VSA'},...
    {'P_g_e_nSA'},{'P_g_e_n_,_s_c_a_l_e_dSA'}];
ytitle = '\DeltaESV_L_V (mL)';
%xtitle = [location,' Negative Work Fraction'];
ymin = -100; ymax = 100;
for i = 1:size(workmetric,1)
    if i == 1 && meas_idx == 2
        xtitle = 'V_fSTNW LBBB';
    elseif i ~= 1 && meas_idx == 2
        xtitle = 'S_fSTNW LBBB';
    elseif i == 1 && meas_idx == 1
        xtitle = 'V_fLVNW LBBB';
    elseif i ~= 1 && meas_idx == 1
        xtitle = 'S_fLVNW LBBB';
    end
       
    fig = plotCorrelationPatientData(deltaESV,workmetric(i,:),titles(i),xtitle,ytitle,ymin,ymax,0,0.6);
    
    %save figures
    % if workmetric == sepnegwork
    %     saveas(fig,['CRT Analysis Plots/esv_corr_SEPnegwork_',figtitle{i},'.png'])
    % elseif workmetric == LVnegwork
    %     saveas(fig,['CRT Analysis Plots/esv_corr_LVnegwork_',figtitle{i},'.png'])
    % end
         [r(i), p(i)] = corr(deltaESV',workmetric(i,:)');
end

p_fisher = FisherRtoZtransformation(r(1),r(2:end),length(deltaESV));


%% LV Negative Work large boxplot
% LVnegwork = [lvNegFractionPSM;lvNegFractionMWCT;LVnegworkMWCT;lvNegFractionMWCTLP;LVnegworkMWCTLP;LVnegworkMWCTLPTV;...
%     LVnegworkMWCT_genericLVP;LVnegworkMWCT_genericLVP_patES];
LVnegwork = [lvNegFractionPSM;lvNegFractionMWCT;lvNegFractionMWCTLP;lvNegFractionMWCTLPTV;...
    lvNegFractionMWCT_genericLVP;lvNegFractionMWCT_genericLVP_patES];
resultPSM = LVnegwork(1,:);%cov_workPSM;%lvNegFractionPSM;
resultCT = LVnegwork(2:end,:);
ymin = -0.1; ymax = 0.5;
ytitle = 'LV Negative Work Fraction';
[fig,mean_results] = plotResponderBoxplotAllEsimates(resultPSM,resultCT,ymin,ymax,ytitle);
saveas(fig,[figpath 'CRT Responder Analysis/LVNW.png'])

%% Septal Negative Work large boxplot
% sepnegwork = [sepNegFractionPSM;sepnegworkMWCT;sepnegworkMWCTLP;sepnegworkMWCTLPTV;...
%     sepnegworkMWCT_genericLVP;sepnegworkMWCT_genericLVP_patES];
sepnegwork = [sepNegFractionPSM;sepNegFractionMWCT;sepNegFractionMWCTLP;sepNegFractionMWCTLPTV;...
    sepNegFractionMWCT_genericLVP;sepNegFractionMWCT_genericLVP_patES];
resultPSM = sepnegwork(1,:);%cov_workPSM;%lvNegFractionPSM;
resultCT = sepnegwork(2:end,:);
ymin = -0.1; ymax = 0.9;
ytitle = 'Septal Negative Work Fraction';
[fig,mean_results] = plotResponderBoxplotAllEsimates(resultPSM,resultCT,ymin,ymax,ytitle);
saveas(fig,[figpath 'CRT Responder Analysis/sepNW.png'])
%% COV large boxplot
% cov = [cov_workPSM;cov_workMWCT;cov_workMWCTLP;cov_workMWCTLPTV;...
%     cov_workMWCT_genericLVP;cov_workMWCT_genericLVP_patES];
cov = [cov_workPSM;covMWCT;covMWCTLP;covMWCTLPTV;covMWCTgenericLVP;covMWCTgenericLVP_patES];
resultPSM = cov(1,:);%cov_workPSM;%lvNegFractionPSM;
resultCT = cov(2:end,:);
ymin = 0; ymax = 3;
ytitle = 'COVW';
[fig,mean_results] = plotResponderBoxplotAllEsimates(resultPSM,resultCT,ymin,ymax,ytitle);
saveas(fig,[figpath 'CRT Responder Analysis/COVW.png'])

