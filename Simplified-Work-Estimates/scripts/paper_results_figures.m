clear
clc;

%lines 147-345: Figure S2
% lines 337 - 494: Figure S3
% lines 676 - 773: Figure S1

figpath = '/Users/amandacraine/Documents/ContijochLab/repos/CRT-PLOS-Submission-pre-repo/figures/';
% get data
%more elaborate MWCT Estimates
load all_seg_work_all_pats.mat %this came from calculateTriMW.m
load segMWCT_effrad_allpats.mat %this came from PrincipalCurvatureAnalysis2.m

%simpler estimates with generic LVP curve
load segMWCT_genericLVP_allpats.mat %this came from evaluateLVPestimates.m

%RSCT
load seg_strain_allpats.mat
load all_seg_strain_all_pats.mat

%% get statistics

segments = 1:17;
pats = 1:8;
%reshape the arrays
segStrain_array = [];
segRSCT_array = [];

segWork_array = reshape(segWork_allpats(segments,pats),numel(segWork_allpats(segments,pats)),1);

segMWCT_array = reshape(segMWCT_allpats(segments,pats),numel(segMWCT_allpats(segments,pats)),1);

segMWCTLP_array = reshape(segMWCTLP(segments,pats),numel(segMWCTLP(segments,pats)),1);
segMWCTLPTV_array = reshape(segMWCTLPTV(segments,pats),numel(segMWCTLPTV(segments,pats)),1);

segMWCT_genericLVP_array = reshape(segMWCT_genericLVP(segments,pats),numel(segMWCT_genericLVP(segments,pats)),1);
segMWCT_genericLVP_patES_array = reshape(segMWCT_genericLVP_patES(segments,pats),numel(segMWCT_genericLVP_patES(segments,pats)),1);

% find the correlations and linear fit for each estimate for each patient
for pat = 1:8
    [r_MWCT_pat(pat),p_MWCT_pat(pat)] = corr(segWork_allpats(segments,pat),segMWCT_allpats(segments,pat),'type','Spearman');
    fit_MWCT_pat(:,pat) = polyfit(segWork_allpats(segments,pat),segMWCT_allpats(segments,pat),1);

    % mdl_MWCT_pat = fitlm(segWork_allpats(segments,pat),segMWCT_allpats(segments,pat));
    % r_MWCT_pat(pat) = mdl_MWCT_pat.Rsquared.Ordinary;
    % p_MWCT_pat(pat) = mdl_MWCT_pat.ModelFitVsNullModel.Pvalue;
    % fit_MWCT_pat(:,pat) = [table2array(mdl_MWCT_pat.Coefficients(2,1)) table2array(mdl_MWCT_pat.Coefficients(1,1))];
    % 
    [r_MWCTLP_pat(pat),p_MWCTLP_pat(pat)] = corr(segWork_allpats(segments,pat),segMWCTLP(segments,pat),'type','Spearman');
    fit_MWCTLP_pat(:,pat) = polyfit(segWork_allpats(segments,pat),segMWCTLP(segments,pat),1);
    
    % mdl_MWCTLP_pat = fitlm(segWork_allpats(segments,pat),segMWCTLP(segments,pat));
    % r_MWCTLP_pat(pat) = mdl_MWCTLP_pat.Rsquared.Ordinary;
    % p_MWCTLP_pat(pat) = mdl_MWCTLP_pat.ModelFitVsNullModel.Pvalue;
    % fit_MWCTLP_pat(:,pat) = [table2array(mdl_MWCTLP_pat.Coefficients(2,1)) table2array(mdl_MWCTLP_pat.Coefficients(1,1))];

    [r_MWCTLPTV_pat(pat),p_MWCTLPTV_pat(pat)] = corr(segWork_allpats(segments,pat),segMWCTLPTV(segments,pat),'type','Spearman');
    fit_MWCTLPTV_pat(:,pat) = polyfit(segWork_allpats(segments,pat),segMWCTLPTV(segments,pat),1);

    % mdl_MWCTLPTV_pat = fitlm(segWork_allpats(segments,pat),segMWCTLPTV(segments,pat));
    % r_MWCTLPTV_pat(pat) = mdl_MWCTLPTV_pat.Rsquared.Ordinary;
    % p_MWCTLPTV_pat(pat) = mdl_MWCTLPTV_pat.ModelFitVsNullModel.Pvalue;
    % fit_MWCTLPTV_pat(:,pat) = [table2array(mdl_MWCTLPTV_pat.Coefficients(2,1)) table2array(mdl_MWCTLPTV_pat.Coefficients(1,1))];
    % 
    [r_MWCT_genericLVP_pat(pat),p_MWCT_genericLVP_pat(pat)] = corr(segWork_allpats(segments,pat),segMWCT_genericLVP(segments,pat),'type','Spearman');
    fit_MWCT_genericLVP_pat(:,pat) = polyfit(segWork_allpats(segments,pat),segMWCT_genericLVP(segments,pat),1);

    % mdl_MWCT_genericLVP_pat = fitlm(segWork_allpats(segments,pat),segMWCT_genericLVP(segments,pat));
    % r_MWCT_genericLVP_pat(pat) = mdl_MWCT_genericLVP_pat.Rsquared.Ordinary;
    % p_MWCT_genericLVP_pat(pat) = mdl_MWCT_genericLVP_pat.ModelFitVsNullModel.Pvalue;
    % fit_MWCT_genericLVP_pat(:,pat) = [table2array(mdl_MWCT_genericLVP_pat.Coefficients(2,1)) table2array(mdl_MWCT_genericLVP_pat.Coefficients(1,1))];
    % 
    [r_MWCT_genericLVP_patES_pat(pat),p_MWCT_genericLVP_patES_pat(pat)] = corr(segWork_allpats(segments,pat),segMWCT_genericLVP_patES(segments,pat),'type','Spearman');
    fit_MWCT_genericLVP_patES_pat(:,pat) = polyfit(segWork_allpats(segments,pat),segMWCT_genericLVP_patES(segments,pat),1);   
    % mdl_MWCT_genericLVP_patES_pat = fitlm(segWork_allpats(segments,pat),segMWCT_genericLVP_patES(segments,pat));
    % r_MWCT_genericLVP_patES_pat(pat) = mdl_MWCT_genericLVP_patES_pat.Rsquared.Ordinary;
    % p_MWCT_genericLVP_patES_pat(pat) = mdl_MWCT_genericLVP_patES_pat.ModelFitVsNullModel.Pvalue;
    % fit_MWCT_genericLVP_patES_pat(:,pat) = [table2array(mdl_MWCT_genericLVP_patES_pat.Coefficients(2,1)) table2array(mdl_MWCT_genericLVP_patES_pat.Coefficients(1,1))];
    % 
end

%%%find correlations and linear fit for each estimate for the whole cohort
[r_MWCT,p_MWCT] = corr(segWork_array,segMWCT_array,'type','Spearman');
fit_MWCT = polyfit(segWork_array,segMWCT_array,1);

% mdl_MWCT = fitlm(segWork_array,segMWCT_array);
% r_MWCT = mdl_MWCT.Rsquared.Ordinary;
% p_MWCT = mdl_MWCT.ModelFitVsNullModel.Pvalue;
% fit_MWCT = [table2array(mdl_MWCT.Coefficients(2,1)) table2array(mdl_MWCT.Coefficients(1,1))]; 
% 
[r_MWCTLP,p_MWCTLP] = corr(segWork_array,segMWCTLP_array,'type','Spearman');
fit_MWCTLP = polyfit(segWork_array,segMWCTLP_array,1);
% mdl_MWCTLP = fitlm(segWork_array,segMWCTLP_array);
% r_MWCTLP = mdl_MWCTLP.Rsquared.Ordinary;
% p_MWCTLP = mdl_MWCTLP.ModelFitVsNullModel.Pvalue;
% fit_MWCTLP = [table2array(mdl_MWCTLP.Coefficients(2,1)) table2array(mdl_MWCTLP.Coefficients(1,1))];
 
[r_MWCTLPTV,p_MWCTLPTV] = corr(segWork_array,segMWCTLPTV_array,'type','Spearman');
fit_MWCTLPTV = polyfit(segWork_array,segMWCTLPTV_array,1);
% mdl_MWCTLPTV = fitlm(segWork_array,segMWCTLPTV_array);
% r_MWCTLPTV = mdl_MWCTLPTV.Rsquared.Ordinary;
% p_MWCTLPTV = mdl_MWCTLPTV.ModelFitVsNullModel.Pvalue;
% fit_MWCTLPTV = [table2array(mdl_MWCTLPTV.Coefficients(2,1)) table2array(mdl_MWCTLPTV.Coefficients(1,1))];
 

[r_MWCT_genericLVP,p_MWCT_genericLVP] = corr(segWork_array,segMWCT_genericLVP_array,'type','Spearman');
fit_MWCT_genericLVP = polyfit(segWork_array,segMWCT_genericLVP_array,1);
% mdl_MWCT_generic = fitlm(segWork_array,segMWCT_genericLVP_array);
% r_MWCT_genericLVP = mdl_MWCT_generic.Rsquared.Ordinary;
% p_MWCT_genericLVP = mdl_MWCT_generic.ModelFitVsNullModel.Pvalue;
% fit_MWCT_genericLVP = [table2array(mdl_MWCT_generic.Coefficients(2,1)) table2array(mdl_MWCT_generic.Coefficients(1,1))];
% 

[r_MWCT_genericLVP_patES,p_MWCT_genericLVP_patES] = corr(segWork_array,segMWCT_genericLVP_patES_array,'type','Spearman');
fit_MWCT_genericLVP_patES = polyfit(segWork_array,segMWCT_genericLVP_patES_array,1);
% mdl_MWCT_generic_patES = fitlm(segWork_array,segMWCT_genericLVP_patES_array);
% r_MWCT_genericLVP_patES = mdl_MWCT_generic_patES.Rsquared.Ordinary;
% p_MWCT_genericLVP_patES = mdl_MWCT_generic_patES.ModelFitVsNullModel.Pvalue;
% fit_MWCT_genericLVP_patES = [table2array(mdl_MWCT_generic_patES.Coefficients(2,1)) table2array(mdl_MWCT_generic_patES.Coefficients(1,1))];

% rmse_psa = rms(segMWCT_array - segWork_array);
% rmse_edwssa = rms(segMWCTLP_array - segWork_array);
% rmse_tvwssa = rms(segMWCTLPTV_array - segWork_array);
% rmse_genpsa = rms(segMWCT_genericLVP_array - segWork_array);
% rmse_genpsa_scales = rms(segMWCT_genericLVP_patES_array - segWork_array);
%maxwork = max([segMWCT_array,segWork_array],[],'all');
%% summary statistics 

%work vs MWCT
summaryMWCT_r = [median(r_MWCT_pat) prctile(r_MWCT_pat,25) prctile(r_MWCT_pat,75)];
summaryMWCT_m = [median(fit_MWCT_pat(1,:)) prctile(fit_MWCT_pat(1,:),25) prctile(fit_MWCT_pat(1,:),75)];

%work vs MWCTLP
summaryMWCTLP_r = [median(r_MWCTLP_pat) prctile(r_MWCTLP_pat,25) prctile(r_MWCTLP_pat,75)];
summaryMWCTLP_m = [median(fit_MWCTLP_pat(1,:)) prctile(fit_MWCTLP_pat(1,:),25) prctile(fit_MWCTLP_pat(1,:),75)];

%work vs MWCTLPTV
summaryMWCTLPTV_r = [median(r_MWCTLPTV_pat) prctile(r_MWCTLPTV_pat,25) prctile(r_MWCTLPTV_pat,75)];
summaryMWCTLPTV_m = [median(fit_MWCTLPTV_pat(1,:)) prctile(fit_MWCTLPTV_pat(1,:),25) prctile(fit_MWCTLPTV_pat(1,:),75)];

%work vs MWCT generic LVP
summaryMWCT_genericLVP_r = [median(r_MWCT_genericLVP_pat) prctile(r_MWCT_genericLVP_pat,25) prctile(r_MWCT_genericLVP_pat,75)];
summaryMWCT_genericLVP_m = [median(fit_MWCT_genericLVP_pat(1,:)) prctile(fit_MWCT_genericLVP_pat(1,:),25) prctile(fit_MWCT_genericLVP_pat(1,:),75)];

%work vs MWCT generic LVP w patient ES
summaryMWCT_genericLVP_patES_r = [median(r_MWCT_genericLVP_patES_pat) prctile(r_MWCT_genericLVP_patES_pat,25) prctile(r_MWCT_genericLVP_patES_pat,75)];
summaryMWCT_genericLVP_patES_m = [median(fit_MWCT_genericLVP_patES_pat(1,:)) prctile(fit_MWCT_genericLVP_patES_pat(1,:),25) prctile(fit_MWCT_genericLVP_patES_pat(1,:),75)];



%% Population plots using the shapes from the crt manuscript
%segments = 1:17;
xaxis_label = "PSM Segmental Work (kJ)";
yaxis_label = "Simplified Seg. Work (kPa)";%"CT Segmental Work (kPa)";
figure; set(gcf,'Position',[0 300 1700 1000]);
for pats = 1:8
    for segments = 1:17
    %define marker shape by patient
    if pats == 1
        markershape = "o";
    elseif pats == 2
        markershape = "square";
    elseif pats == 3
        markershape = "^";
    elseif pats == 4
        markershape = "v";
    elseif pats == 5
        markershape = "diamond";
    elseif pats == 6
        markershape = "pentagram";
    elseif pats == 7
        markershape = "hexagram";
    else
        markershape = ">";
    end

    %define marker fill by segment
    %basal segments
    if segments == 1||segments == 2||segments == 3||segments == 4||segments == 5||segments == 6
        markeredge = "k";
        markerfill = "w";
    %mid segments
    elseif segments == 7||segments == 8||segments == 9||segments == 10||segments == 11||segments == 12
        markeredge = [0.5 0.5 0.5];
        markerfill = [0.5 0.5 0.5];
     %apical segments
    else
        markeredge = "k";
        markerfill = "k";
    end
        subplot(2,3,1); hold all; %set(gca,'colororder',turbo(8)); hold all;
    plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
    plot(segWork_allpats(segments,pats),segMWCT_allpats(segments,pats),...
        markershape,'MarkerFaceColor',markerfill,'MarkerEdgeColor',markeredge)
    if pats == 1 && segments == 1
        text(3.7,-2,['R = ',num2str(r_MWCT,'%.2f')],'HorizontalAlignment','right')
        if p_MWCT > 0.01
            text(3.7,-3,['p = ',num2str(p_MWCT,'%.2f')],'HorizontalAlignment','right')
        else
            text(3.7,-3,'p < 0.01','HorizontalAlignment','right')
        end
        text(3.7,-4,['y = ',num2str(fit_MWCT(1),'%.2f'),'x + ',num2str(fit_MWCT(2),'%.2f')],'HorizontalAlignment','right')
        % else
        %     text(3.7,-2,['r = ',num2str(r_MWCT_pat(pats),'%.2f')],'HorizontalAlignment','right')
        %     if p_MWCT_pat(pats) > 0.01
        %         text(3.7,-3,['p = ',num2str(p_MWCT_pat(pats),'%.2f')],'HorizontalAlignment','right')
        %     else
        %         text(3.7,-3,'p < 0.01','HorizontalAlignment','right')
        %     end
        %     text(3.7,-4,['y = ',num2str(fit_MWCT_pat(1,pats),'%.2f'),'x + ',num2str(fit_MWCT_pat(2,pats),'%.2f')],'HorizontalAlignment','right')
    end
    axis([-0.5 4 -4.5 4],'square')
    xlabel(xaxis_label)
    ylabel(yaxis_label)
    %title('PSA_L_H_C') %Pressure-strain area, left heart cath
    title('P_L_H_CSA')

    subplot(2,3,2); hold all;% set(gca,'colororder',turbo(8)); hold all;
    plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
    plot(segWork_allpats(segments,pats),segMWCTLP(segments,pats),...
        markershape,'MarkerFaceColor',markerfill,'MarkerEdgeColor',markeredge)
    if pats == 1 && segments == 1 %numel(pats) > 1
        text(3.7,-2,['R = ',num2str(r_MWCTLP,'%.2f')],'HorizontalAlignment','right')
        if p_MWCTLP > 0.01
            text(3.7,-3,['p = ',num2str(p_MWCTLP,'%.2f')],'HorizontalAlignment','right')
        else
            text(3.7,-3,'p < 0.01','HorizontalAlignment','right')
        end
        text(3.7,-4,['y = ',num2str(fit_MWCTLP(1),'%.2f'),'x + ',num2str(fit_MWCTLP(2),'%.2f')],'HorizontalAlignment','right')
    % else
    %     text(3.7,-2,['r = ',num2str(r_MWCTLP_pat(pats),'%.2f')],'HorizontalAlignment','right')
    %     if p_MWCTLP_pat(pats) > 0.01
    %         text(3.7,-3,['p = ',num2str(p_MWCTLP_pat(pats),'%.2f')],'HorizontalAlignment','right')
    %     else
    %         text(3.7,-3,'p < 0.01','HorizontalAlignment','right')
    %     end
    %     text(3.7,-4,['y = ',num2str(fit_MWCTLP_pat(1,pats),'%.2f'),'x + ',num2str(fit_MWCTLP_pat(2,pats),'%.2f')],'HorizontalAlignment','right')
    end
    axis([-0.5 4 -4.5 4],'square')
    xlabel(xaxis_label)
    ylabel(yaxis_label)
    %title('WSSA_E_D') %Wall Stress-strain area, end-diastole
    title('WS_E_DSA')

    subplot(2,3,3); hold all; %set(gca,'colororder',turbo(8)); hold all;
    plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
    plot(segWork_allpats(segments,pats),segMWCTLPTV(segments,pats),...
        markershape,'MarkerFaceColor',markerfill,'MarkerEdgeColor',markeredge)
    if pats == 1 && segments == 1 %numel(pats) > 1
        text(3.7,-2,['R = ',num2str(r_MWCTLPTV,'%.2f')],'HorizontalAlignment','right')
        if p_MWCTLPTV > 0.01
            text(3.7,-3,['p = ',num2str(p_MWCTLPTV,'%.2f')],'HorizontalAlignment','right')
        else
            text(3.7,-3,'p < 0.01','HorizontalAlignment','right')
        end
        text(3.7,-4,['y = ',num2str(fit_MWCTLPTV(1),'%.2f'),'x + ',num2str(fit_MWCTLPTV(2),'%.2f')],'HorizontalAlignment','right')
    % else
    %     text(3.7,-2,['r = ',num2str(r_MWCTLPTV_pat(pats),'%.2f')],'HorizontalAlignment','right')
    %     if p_MWCTLPTV_pat(pats) > 0.01
    %         text(3.7,-3,['p = ',num2str(p_MWCTLPTV_pat(pats),'%.2f')],'HorizontalAlignment','right')
    %     else
    %         text(3.7,-3,'p < 0.01','HorizontalAlignment','right')
    %     end
    %     text(3.7,-4,['y = ',num2str(fit_MWCTLPTV_pat(1,pats),'%.2f'),'x + ',num2str(fit_MWCTLPTV_pat(2,pats),'%.2f')],'HorizontalAlignment','right')
    end
    axis([-0.5 4 -4.5 4],'square')
    xlabel(xaxis_label)
    ylabel(yaxis_label)
    %title('WSSA_T_V') %Wall stress-strain area, time-varying
    title('WS_T_VSA')

    subplot(2,3,4); hold all; %set(gca,'colororder',turbo(8)); hold all;
    plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
    plot(segWork_allpats(segments,pats),segMWCT_genericLVP(segments,pats),...
        markershape,'MarkerFaceColor',markerfill,'MarkerEdgeColor',markeredge)
    if pats == 1 && segments == 1%numel(pats) > 1
        text(3.7,-2,['R = ',num2str(r_MWCT_genericLVP,'%.2f')],'HorizontalAlignment','right')
        if p_MWCT_genericLVP > 0.01
            text(3.7,-3,['p = ',num2str(p_MWCT_genericLVP,'%.2f')],'HorizontalAlignment','right')
        else
            text(3.7,-3,'p < 0.01','HorizontalAlignment','right')
        end
        text(3.7,-4,['y = ',num2str(fit_MWCT_genericLVP(1),'%.2f'),'x + ',num2str(fit_MWCT_genericLVP(2),'%.2f')],'HorizontalAlignment','right')
    % else
    %     text(3.7,-2,['r = ',num2str(r_MWCT_genericLVP_pat(pats),'%.2f')],'HorizontalAlignment','right')
    %     if p_MWCT_genericLVP_pat(pats) > 0.01
    %         text(3.7,-3,['p = ',num2str(p_MWCT_genericLVP_pat(pats),'%.2f')],'HorizontalAlignment','right')
    %     else
    %         text(3.7,-3,'p < 0.01','HorizontalAlignment','right')
    %     end
    %     text(3.7,-4,['y = ',num2str(fit_MWCT_genericLVP_pat(1,pats),'%.2f'),'x + ',num2str(fit_MWCT_genericLVP_pat(2,pats),'%.2f')],'HorizontalAlignment','right')
    end
    axis([-0.5 4 -4.5 4],'square')
    xlabel(xaxis_label)
    ylabel(yaxis_label)
    %title('PSA_g_e_n') %pressure-strain area, generic pressure
    title('P_g_e_nSA')

    subplot(2,3,5); hold all; set(gca,'colororder',turbo(8)); hold all;
    plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
    plot(segWork_allpats(segments,pats),segMWCT_genericLVP_patES(segments,pats),...
        markershape,'MarkerFaceColor',markerfill,'MarkerEdgeColor',markeredge)
    if pats == 1 && segments == 1%numel(pats) > 1
        text(3.7,-2,['R = ',num2str(r_MWCT_genericLVP_patES,'%.2f')],'HorizontalAlignment','right')
        if p_MWCT_genericLVP_patES > 0.01
            text(3.7,-3,['p = ',num2str(p_MWCT_genericLVP_patES,'%.2f')],'HorizontalAlignment','right')
        else
            text(3.7,-3,'p < 0.01','HorizontalAlignment','right')
        end
        text(3.7,-4,['y = ',num2str(fit_MWCT_genericLVP_patES(1),'%.2f'),'x + ',num2str(fit_MWCT_genericLVP_patES(2),'%.2f')],'HorizontalAlignment','right')
    % else
    %     text(3.7,-2,['r = ',num2str(r_MWCT_genericLVP_patES_pat(pats),'%.2f')],'HorizontalAlignment','right')
    %     if p_MWCT_genericLVP_patES_pat(pats) > 0.01
    %         text(3.7,-3,['p = ',num2str(p_MWCT_genericLVP_patES_pat(pats),'%.2f')],'HorizontalAlignment','right')
    %     else
    %         text(3.7,-3,'p < 0.01','HorizontalAlignment','right')
    %     end
    %     text(3.7,-4,['y = ',num2str(fit_MWCT_genericLVP_patES_pat(1,pats),'%.2f'),'x + ',num2str(fit_MWCT_genericLVP_patES_pat(2,pats),'%.2f')],'HorizontalAlignment','right')
    end
    axis([-0.5 4 -4.5 4],'square')
    xlabel(xaxis_label)
    ylabel(yaxis_label)
    %title('PSA_g_e_n_,_s_c_a_l_e_d') %pressure-strain area, generic pressure scaled to peak pressure
    title('P_g_e_n_,_s_c_a_l_e_dSA')
    % if numel(pats) > 1
    %     legend('','Pat 1','Pat 2','Pat 3','Pat 4','Pat 5','Pat 6','Pat 7','Pat 8','position',[0.95 0.69 0 0])
    % end

    
    end
    
end

%sgtitle('Agreement between Segmental PSM Work and Segmental LV Work Estimates')

set(findall(gcf,'-property','FontSize'),'FontSize',22)
set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',10)
set(findall(gcf,'-property','Box'),'Box','on')

saveas(gcf,[figpath 'Supplement Figures/S2_seg_work_pop.png'])

%% Population Plots, one patient
%segments = 1:17;
xaxis_label = "PSM Segmental Work (kJ)";
yaxis_label = "Simplified Seg. Work (kPa)";%"CT Segmental Work (kPa)";
pats = 5;
if pats == 1
    markershape = "o";
elseif pats == 2
    markershape = "square";
elseif pats == 3
    markershape = "^";
elseif pats == 4
    markershape = "v";
elseif pats == 5
    markershape = "diamond";
elseif pats == 6
    markershape = "pentagram";
elseif pats == 7
    markershape = "hexagram";
else
    markershape = ">";
end

figure; set(gcf,'Position',[0 300 1700 1000]);
for segments = 1:17
    %define marker fill by segment
    %basal segments
    if segments == 1||segments == 2||segments == 3||segments == 4||segments == 5||segments == 6
        markeredge = "k";
        markerfill = "w";
        %mid segments
    elseif segments == 7||segments == 8||segments == 9||segments == 10||segments == 11||segments == 12
        markeredge = [0.5 0.5 0.5];
        markerfill = [0.5 0.5 0.5];
        %apical segments
    else
        markeredge = "k";
        markerfill = "k";
    end

    %figure; set(gcf,'Position',[0 300 1700 1000]);
    subplot(2,3,1); hold all; %set(gca,'colororder',turbo(8)); hold all;
    plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
    plot(segWork_allpats(segments,pats),segMWCT_allpats(segments,pats),...
        markershape,'MarkerFaceColor',markerfill,'MarkerEdgeColor',markeredge)
    if segments == 17
        text(2.3,-0.2,['R = ',num2str(r_MWCT_pat(pats),'%.2f')],'HorizontalAlignment','right')
        if p_MWCT_pat(pats) > 0.01
            text(2.3,-0.5,['p = ',num2str(p_MWCT_pat(pats),'%.2f')],'HorizontalAlignment','right')
        else
            text(2.3,-0.5,'p < 0.01','HorizontalAlignment','right')
        end
        text(2.3,-0.8,['y = ',num2str(fit_MWCT_pat(1,pats),'%.2f'),'x + ',num2str(fit_MWCT_pat(2,pats),'%.2f')],'HorizontalAlignment','right')
    end
    axis([-1 2.5 -1 2.5],'square')
    xlabel(xaxis_label)
    ylabel(yaxis_label)
    title('P_L_H_CSA')

    subplot(2,3,2); hold all; %set(gca,'colororder',turbo(8)); hold all;
    plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
    plot(segWork_allpats(segments,pats),segMWCTLP(segments,pats),...
        markershape,'MarkerFaceColor',markerfill,'MarkerEdgeColor',markeredge)
    if segments == 17
        text(2.3,-0.2,['R = ',num2str(r_MWCTLP_pat(pats),'%.2f')],'HorizontalAlignment','right')
    
        if p_MWCTLP_pat(pats) > 0.01
            text(2.3,-0.5,['p = ',num2str(p_MWCTLP_pat(pats),'%.2f')],'HorizontalAlignment','right')
        else
            text(2.3,-0.5,'p < 0.01','HorizontalAlignment','right')
        end
        text(2.3,-0.8,['y = ',num2str(fit_MWCTLP_pat(1,pats),'%.2f'),'x + ',num2str(fit_MWCTLP_pat(2,pats),'%.2f')],'HorizontalAlignment','right')
    end
    axis([-1 2.5 -1 2.5],'square')
    xlabel(xaxis_label)
    ylabel(yaxis_label)
    title('WS_E_DSA')

    subplot(2,3,3); hold all; % set(gca,'colororder',turbo(8)); hold all;
    plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
    plot(segWork_allpats(segments,pats),segMWCTLPTV(segments,pats),...
        markershape,'MarkerFaceColor',markerfill,'MarkerEdgeColor',markeredge)
    if segments == 17
        text(2.3,-0.2,['R = ',num2str(r_MWCTLPTV_pat(pats),'%.2f')],'HorizontalAlignment','right')
    
        if p_MWCTLPTV_pat(pats) > 0.01
            text(2.3,-0.5,['p = ',num2str(p_MWCTLPTV_pat(pats),'%.2f')],'HorizontalAlignment','right')
        else
            text(2.3,-0.5,'p < 0.01','HorizontalAlignment','right')
        end
        text(2.3,-0.8,['y = ',num2str(fit_MWCTLPTV_pat(1,pats),'%.2f'),'x + ',num2str(fit_MWCTLPTV_pat(2,pats),'%.2f')],'HorizontalAlignment','right')
    end
    axis([-1 2.5 -1 2.5],'square')
    xlabel(xaxis_label)
    ylabel(yaxis_label)
    title('WS_T_VSA')

    subplot(2,3,4); hold all; %set(gca,'colororder',turbo(8)); hold all;
    plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
    plot(segWork_allpats(segments,pats),segMWCT_genericLVP(segments,pats),...
        markershape,'MarkerFaceColor',markerfill,'MarkerEdgeColor',markeredge)
    if segments == 17
        text(2.3,-0.2,['R = ',num2str(r_MWCT_genericLVP_pat(pats),'%.2f')],'HorizontalAlignment','right')
    
        if p_MWCT_genericLVP_pat(pats) > 0.01
            text(2.3,-0.5,['p = ',num2str(p_MWCT_genericLVP_pat(pats),'%.2f')],'HorizontalAlignment','right')
        else
            text(2.3,-0.5,'p < 0.01','HorizontalAlignment','right')
        end
        if fit_MWCT_genericLVP_pat > 0
            text(2.3,-0.8,['y = ',num2str(fit_MWCT_genericLVP_pat(1,pats),'%.2f'),'x + ',num2str(fit_MWCT_genericLVP_pat(2,pats),'%.2f')],'HorizontalAlignment','right')
        else
            text(2.3,-0.8,['y = ',num2str(fit_MWCT_genericLVP_pat(1,pats),'%.2f'),'x - ',num2str(abs(fit_MWCT_genericLVP_pat(2,pats)),'%.2f')],'HorizontalAlignment','right')
        end
    end
    axis([-1 2.5 -1 2.5],'square')
    xlabel(xaxis_label)
    ylabel(yaxis_label)
    title('P_g_e_nSA')

    subplot(2,3,5); hold all; %set(gca,'colororder',turbo(8)); hold all;
    plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
    plot(segWork_allpats(segments,pats),segMWCT_genericLVP_patES(segments,pats),...
        markershape,'MarkerFaceColor',markerfill,'MarkerEdgeColor',markeredge)
    if segments == 17
        text(2.3,-0.2,['R = ',num2str(r_MWCT_genericLVP_patES_pat(pats),'%.2f')],'HorizontalAlignment','right')
    
        if p_MWCT_genericLVP_patES_pat(pats) > 0.01
            text(2.3,-0.5,['p = ',num2str(p_MWCT_genericLVP_patES_pat(pats),'%.2f')],'HorizontalAlignment','right')
        else
            text(2.3,-0.5,'p < 0.01','HorizontalAlignment','right')
        end
        if fit_MWCT_genericLVP_patES_pat(2,pats) > 0
            text(2.3,-0.8,['y = ',num2str(fit_MWCT_genericLVP_patES_pat(1,pats),'%.2f'),'x + ',num2str(fit_MWCT_genericLVP_patES_pat(2,pats),'%.2f')],'HorizontalAlignment','right')
        else
            text(2.3,-0.8,['y = ',num2str(fit_MWCT_genericLVP_patES_pat(1,pats),'%.2f'),'x - ',num2str(abs(fit_MWCT_genericLVP_patES_pat(2,pats)),'%.2f')],'HorizontalAlignment','right')
        end
    end
    axis([-1 2.5 -1 2.5],'square')
    xlabel(xaxis_label)
    ylabel(yaxis_label)
    title('P_g_e_n_,_s_c_a_l_e_dSA')

end

if numel(pats) > 1
    legend('','Pat 1','Pat 2','Pat 3','Pat 4','Pat 5','Pat 6','Pat 7','Pat 8','position',[0.95 0.69 0 0])
end


%sgtitle('Agreement between Segmental PSM Work and Segmental LV Work Estimates, Patient 5')

set(findall(gcf,'-property','FontSize'),'FontSize',22)
set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',10)
set(findall(gcf,'-property','Box'),'Box','on')

saveas(gcf,[figpath 'Supplement Figures/S3_seg_work_pat.png'])

%% Population plots using the shapes from the crt manuscript, mean LV Work
segments = 1:17;
figure; set(gcf,'Position',[0 300 1700 1000]);
for pats = 1:8
    %figure; set(gcf,'Position',[0 300 1700 1000]);
    if pats == 1 || pats == 3 || pats == 5 || pats == 6
        markerfill = 'w';
    elseif pats == 2 || pats == 4 || pats == 8
        markerfill = [0.5 0.5 0.5];
    else
        markerfill = 'k';
    end
    if pats == 2 || pats == 4 || pats == 5
        markershape = "diamond";
    elseif pats == 1 || pats == 6 || pats == 7
        markershape = "o";
    else
        markershape = "^";
    end
    if pats == 2 || pats == 4 || pats == 8
        markeredge = [0.5 0.5 0.5];
    else
        markeredge = 'k';
    end

    subplot(2,3,1); hold all; %set(gca,'colororder',turbo(8)); hold all;
    plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
    plot(mean(segWork_allpats(segments,pats)),mean(segMWCT_allpats(segments,pats)),...
        markershape,'MarkerFaceColor',markerfill,'MarkerEdgeColor',markeredge)
    [r_MWCT_mean,p_MWCT_mean] = corr(mean(segWork_allpats)',mean(segMWCT_allpats)','type','Spearman');
    fit_MWCT_mean = polyfit(mean(segWork_allpats)',mean(segMWCT_allpats)',1);
    if pats == 8
        text(1.8,0.5,['r = ',num2str(r_MWCT_mean,'%.2f')],'HorizontalAlignment','right')
        if p_MWCT_mean > 0.01
            text(1.8,0.3,['p = ',num2str(p_MWCT_mean,'%.2f')],'HorizontalAlignment','right')
        else
            text(1.8,0.3,'p < 0.01','HorizontalAlignment','right')
        end
        text(1.8,0.1,['y = ',num2str(fit_MWCT_mean(1),'%.2f'),'x + ',num2str(fit_MWCT_mean(2),'%.2f')],'HorizontalAlignment','right')
        % else
        %     text(3.7,-2,['r = ',num2str(r_MWCT_pat(pats),'%.2f')],'HorizontalAlignment','right')
        %     if p_MWCT_pat(pats) > 0.01
        %         text(3.7,-3,['p = ',num2str(p_MWCT_pat(pats),'%.2f')],'HorizontalAlignment','right')
        %     else
        %         text(3.7,-3,'p < 0.01','HorizontalAlignment','right')
        %     end
        %     text(3.7,-4,['y = ',num2str(fit_MWCT_pat(1,pats),'%.2f'),'x + ',num2str(fit_MWCT_pat(2,pats),'%.2f')],'HorizontalAlignment','right')
    end
    axis([0 2 0 2],'square')
    xlabel('LV Segmental Work via PSM (kJ)')
    ylabel('LV Segmental Work via CT (kPa)')
    title('Catheter LV Pressure')

    subplot(2,3,2); hold all;% set(gca,'colororder',turbo(8)); hold all;
    plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
    plot(mean(segWork_allpats(segments,pats)),mean(segMWCTLP(segments,pats)),...
        markershape,'MarkerFaceColor',markerfill,'MarkerEdgeColor',markeredge)
    [r_MWCTLP_mean,p_MWCTLP_mean] = corr(mean(segWork_allpats)',mean(segMWCTLP)','type','Spearman');
    fit_MWCTLP_mean = polyfit(mean(segWork_allpats)',mean(segMWCTLP)',1);
    if pats == 8
        text(1.8,0.5,['r = ',num2str(r_MWCTLP_mean,'%.2f')],'HorizontalAlignment','right')
        if p_MWCTLP_mean > 0.01
            text(1.8,0.3,['p = ',num2str(p_MWCTLP_mean,'%.2f')],'HorizontalAlignment','right')
        else
            text(1.8,0.3,'p < 0.01','HorizontalAlignment','right')
        end
        text(1.8,0.1,['y = ',num2str(fit_MWCTLP_mean(1),'%.2f'),'x + ',num2str(fit_MWCTLP_mean(2),'%.2f')],'HorizontalAlignment','right')
        
    % else
    %     text(3.7,-2,['r = ',num2str(r_MWCTLP_pat(pats),'%.2f')],'HorizontalAlignment','right')
    %     if p_MWCTLP_pat(pats) > 0.01
    %         text(3.7,-3,['p = ',num2str(p_MWCTLP_pat(pats),'%.2f')],'HorizontalAlignment','right')
    %     else
    %         text(3.7,-3,'p < 0.01','HorizontalAlignment','right')
    %     end
    %     text(3.7,-4,['y = ',num2str(fit_MWCTLP_pat(1,pats),'%.2f'),'x + ',num2str(fit_MWCTLP_pat(2,pats),'%.2f')],'HorizontalAlignment','right')
    end
    axis([0 2 0 2],'square')
    xlabel('LV Segmental Work via PSM (kJ)')
    ylabel('LV Segmental Work via CT (kPa)')
    title('Static Laplace Wall Stress')

    subplot(2,3,3); hold all; %set(gca,'colororder',turbo(8)); hold all;
    plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
    plot(mean(segWork_allpats(segments,pats)),mean(segMWCTLPTV(segments,pats)),...
        markershape,'MarkerFaceColor',markerfill,'MarkerEdgeColor',markeredge)
    [r_MWCTLPTV_mean,p_MWCTLPTV_mean] = corr(mean(segWork_allpats)',mean(segMWCTLPTV)','type','Spearman');
    fit_MWCTLPTV_mean = polyfit(mean(segWork_allpats)',mean(segMWCTLPTV)',1);
    if pats == 8
        text(1.8,0.5,['r = ',num2str(r_MWCTLPTV_mean,'%.2f')],'HorizontalAlignment','right')
        if p_MWCTLPTV_mean > 0.01
            text(1.8,0.3,['p = ',num2str(p_MWCTLPTV_mean,'%.2f')],'HorizontalAlignment','right')
        else
            text(1.8,0.3,'p < 0.01','HorizontalAlignment','right')
        end
        text(1.8,0.1,['y = ',num2str(fit_MWCTLPTV_mean(1),'%.2f'),'x + ',num2str(fit_MWCTLPTV_mean(2),'%.2f')],'HorizontalAlignment','right')
        
    % else
    %     text(3.7,-2,['r = ',num2str(r_MWCTLPTV_pat(pats),'%.2f')],'HorizontalAlignment','right')
    %     if p_MWCTLPTV_pat(pats) > 0.01
    %         text(3.7,-3,['p = ',num2str(p_MWCTLPTV_pat(pats),'%.2f')],'HorizontalAlignment','right')
    %     else
    %         text(3.7,-3,'p < 0.01','HorizontalAlignment','right')
    %     end
    %     text(3.7,-4,['y = ',num2str(fit_MWCTLPTV_pat(1,pats),'%.2f'),'x + ',num2str(fit_MWCTLPTV_pat(2,pats),'%.2f')],'HorizontalAlignment','right')
    end
    axis([0 2 0 2],'square')
    xlabel('LV Segmental Work via PSM (kJ)')
    ylabel('LV Segmental Work via CT (kPa)')
    title([{'Time-Varying'};{' Laplace Wall Stress'}])

    subplot(2,3,4); hold all; %set(gca,'colororder',turbo(8)); hold all;
    plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
    plot(mean(segWork_allpats(segments,pats)),mean(segMWCT_genericLVP(segments,pats)),...
        markershape,'MarkerFaceColor',markerfill,'MarkerEdgeColor',markeredge)
    [r_MWCT_genericLVP_mean,p_MWCT_genericLVP_mean] = corr(mean(segWork_allpats)',mean(segMWCT_genericLVP)','type','Spearman');
    fit_MWCT_genericLVP_mean = polyfit(mean(segWork_allpats)',mean(segMWCT_genericLVP)',1);
    if pats == 8
        text(1.8,0.5,['r = ',num2str(r_MWCT_genericLVP_mean,'%.2f')],'HorizontalAlignment','right')
        if p_MWCT_genericLVP_mean > 0.01
            text(1.8,0.3,['p = ',num2str(p_MWCT_genericLVP_mean,'%.2f')],'HorizontalAlignment','right')
        else
            text(1.8,0.3,'p < 0.01','HorizontalAlignment','right')
        end
        text(1.8,0.1,['y = ',num2str(fit_MWCT_genericLVP_mean(1),'%.2f'),'x + ',num2str(fit_MWCT_genericLVP_mean(2),'%.2f')],'HorizontalAlignment','right')
    % else
    %     text(3.7,-2,['r = ',num2str(r_MWCT_genericLVP_pat(pats),'%.2f')],'HorizontalAlignment','right')
    %     if p_MWCT_genericLVP_pat(pats) > 0.01
    %         text(3.7,-3,['p = ',num2str(p_MWCT_genericLVP_pat(pats),'%.2f')],'HorizontalAlignment','right')
    %     else
    %         text(3.7,-3,'p < 0.01','HorizontalAlignment','right')
    %     end
    %     text(3.7,-4,['y = ',num2str(fit_MWCT_genericLVP_pat(1,pats),'%.2f'),'x + ',num2str(fit_MWCT_genericLVP_pat(2,pats),'%.2f')],'HorizontalAlignment','right')
    end
    axis([0 2 0 2],'square')
    xlabel('LV Segmental Work via PSM (kJ)')
    ylabel('LV Segmental Work via CT (kPa)')
    title('Generic LV Pressure')

    subplot(2,3,5); hold all; set(gca,'colororder',turbo(8)); hold all;
    plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
    plot(mean(segWork_allpats(segments,pats)),mean(segMWCT_genericLVP_patES(segments,pats)),...
        markershape,'MarkerFaceColor',markerfill,'MarkerEdgeColor',markeredge)
    [r_MWCT_genericLVP_patES_mean,p_MWCT_genericLVP_patES_mean] = corr(mean(segWork_allpats)',mean(segMWCT_genericLVP_patES)','type','Spearman');
    fit_MWCT_genericLVP_patES_mean = polyfit(mean(segWork_allpats)',mean(segMWCT_genericLVP_patES)',1);
    if pats == 8
        text(1.8,0.5,['r = ',num2str(r_MWCT_genericLVP_patES_mean,'%.2f')],'HorizontalAlignment','right')
        if p_MWCT_genericLVP_patES_mean > 0.01
            text(1.8,0.3,['p = ',num2str(p_MWCT_genericLVP_patES_mean,'%.2f')],'HorizontalAlignment','right')
        else
            text(1.8,0.3,'p < 0.01','HorizontalAlignment','right')
        end
        text(1.8,0.1,['y = ',num2str(fit_MWCT_genericLVP_patES_mean(1),'%.2f'),'x + ',num2str(fit_MWCT_genericLVP_patES_mean(2),'%.2f')],'HorizontalAlignment','right')
    % else
    %     text(3.7,-2,['r = ',num2str(r_MWCT_genericLVP_patES_pat(pats),'%.2f')],'HorizontalAlignment','right')
    %     if p_MWCT_genericLVP_patES_pat(pats) > 0.01
    %         text(3.7,-3,['p = ',num2str(p_MWCT_genericLVP_patES_pat(pats),'%.2f')],'HorizontalAlignment','right')
    %     else
    %         text(3.7,-3,'p < 0.01','HorizontalAlignment','right')
    %     end
    %     text(3.7,-4,['y = ',num2str(fit_MWCT_genericLVP_patES_pat(1,pats),'%.2f'),'x + ',num2str(fit_MWCT_genericLVP_patES_pat(2,pats),'%.2f')],'HorizontalAlignment','right')
    end
    axis([0 2 0 2],'square')
    xlabel('LV Segmental Work via PSM (kJ)')
    ylabel('LV Segmental Work via CT (kPa)')
    title([{'Generic LV Pressure with'};{'Patient Peak Pressure'}])
    % if numel(pats) > 1
    %     legend('','Pat 1','Pat 2','Pat 3','Pat 4','Pat 5','Pat 6','Pat 7','Pat 8','position',[0.95 0.69 0 0])
    % end

end

sgtitle('Agreement between Segmental PSM Work and Segmental LV Work Estimates')

set(findall(gcf,'-property','FontSize'),'FontSize',22)
set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',10)
set(findall(gcf,'-property','Box'),'Box','on')


%% Population plots using the shapes from the crt manuscript, strain
%segments = 1:17;
for pats = 1:8
    peakstrain(:,pats) = min(seg_strain_allpats{pats},[],2);
    peakRSCT(:,pats) = min(segRSCT_allpats{pats},[],2);

    [r_rsct_pats(pats),p_rsct_pats(pats)] = corr(peakstrain(:,pats),peakRSCT(:,pats),'type','Spearman');
    rsct_fit_pats(:,pats) = polyfit(peakstrain(:,pats),peakRSCT(:,pats),1);
    % mdl = fitlm(peakstrain(:,pats),peakRSCT(:,pats));
    % 
    % r_rsct_pats(pats) = mdl.Rsquared.Ordinary;
    % p_rsct_pats(pats) = mdl.ModelFitVsNullModel.Pvalue;
    % rsct_fit_pats(:,pats) = [table2array(mdl.Coefficients(2,1)) table2array(mdl.Coefficients(1,1))];
    % 

end


peakstrain_array = reshape(peakstrain,numel(peakstrain),1);
peakRSCT_array = reshape(peakRSCT,numel(peakRSCT),1);

[r_rsct,p_rsct] = corr(peakstrain_array,peakRSCT_array,'type','Spearman');
% mdl_rsct = fitlm(peakstrain_array,peakRSCT_array);
% r_rsct = mdl_rsct.Rsquared.Ordinary;
% p_rsct = mdl_rsct.ModelFitVsNullModel.Pvalue;
% fit_rsct = [table2array(mdl_rsct.Coefficients(2,1)) table2array(mdl_rsct.Coefficients(1,1))];
fit_rsct = polyfit(peakstrain_array,peakRSCT_array,1);

figure; set(gcf,'Position',[300 300 600 600]);
for pats = 1:8
    for segments = 1:17
    %define marker shape by patient
    if pats == 1
        markershape = "o";
    elseif pats == 2
        markershape = "square";
    elseif pats == 3
        markershape = "^";
    elseif pats == 4
        markershape = "v";
    elseif pats == 5
        markershape = "diamond";
    elseif pats == 6
        markershape = "pentagram";
    elseif pats == 7
        markershape = "hexagram";
    else
        markershape = ">";
    end

    %define marker fill by segment
    %basal segments
    if segments == 1||segments == 2||segments == 3||segments == 4||segments == 5||segments == 6
        markeredge = "k";
        markerfill = "w";
    %mid segments
    elseif segments == 7||segments == 8||segments == 9||segments == 10||segments == 11||segments == 12
        markeredge = [0.5 0.5 0.5];
        markerfill = [0.5 0.5 0.5];
     %apical segments
    else
        markeredge = "k";
        markerfill = "k";
    end
     hold all;  %set(gca,'colororder',turbo(8)); hold all;
    plot([-0.2 0.02],[ -0.2 0.02],'--','Color',[0.9 0.9 0.9])
    plot(peakstrain(segments,pats),peakRSCT(segments,pats),...
        markershape,'MarkerFaceColor',markerfill,'MarkerEdgeColor',markeredge)
    if pats == 1 && segments == 1
        text(0,-0.15,['R = ',num2str(r_rsct,'%.2f')],'HorizontalAlignment','right')
        if p_rsct > 0.01
            text(0,-0.165,['p = ',num2str(p_rsct,'%.2f')],'HorizontalAlignment','right')
        else
            text(0,-0.165,'p < 0.01','HorizontalAlignment','right')
        end
        if fit_rsct(2) >= 0
            text(0,-0.18,['y = ',num2str(fit_rsct(1),'%.2f'),'x + ',num2str(fit_rsct(2),'%.2f')],'HorizontalAlignment','right')
        else
            text(0,-0.18,['y = ',num2str(fit_rsct(1),'%.2f'),'x - ',num2str(abs(fit_rsct(2)),'%.2f')],'HorizontalAlignment','right')
        end
    end
    axis([-0.2 0.02 -0.2 0.02],'square')
    xlabel('PSM Segmental Peak Strain')
    %ylabel('CT Segmental Peak Strain')
    ylabel('Simplified Seg. Peak Strain')
    title('Segmental LV Peak Strain') 
    %title('P_L_H_CSA')
    end
end



%sgtitle('Agreement between Segmental PSM Work and Segmental LV Work Estimates')

set(findall(gcf,'-property','FontSize'),'FontSize',22)
set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',10)
set(findall(gcf,'-property','Box'),'Box','on')

saveas(gcf,[figpath 'Supplement Figures/S1_seg_strain_pop.png'])
%% Population Plots, one patient
%segments = 1:17;
for pats = 1:8
    peakstrain(:,pats) = min(seg_strain_allpats{pats},[],2);
    peakRSCT(:,pats) = min(segRSCT_allpats{pats},[],2);

    [r_rsct_pats(pats),p_rsct_pats(pats)] = corr(peakstrain(:,pats),peakRSCT(:,pats),'type','Spearman');
    rsct_fit_pats(:,pats) = polyfit(peakstrain(:,pats),peakRSCT(:,pats),1);
end
peakstrain_array = reshape(peakstrain,numel(peakstrain),1);
peakRSCT_array = reshape(peakRSCT,numel(peakRSCT),1);

[r_rsct,p_rsct] = corr(peakstrain_array,peakRSCT_array,'type','Spearman');
fit_rsct = polyfit(peakstrain_array,peakRSCT_array,1);

xaxis_label = "PSM Segmental Work (kJ)";
yaxis_label = "CT Segmental Work (kPa)";
pats = 8;
if pats == 1
    markershape = "o";
elseif pats == 2
    markershape = "square";
elseif pats == 3
    markershape = "^";
elseif pats == 4
    markershape = "v";
elseif pats == 5
    markershape = "diamond";
elseif pats == 6
    markershape = "pentagram";
elseif pats == 7
    markershape = "hexagram";
else
    markershape = ">";
end

figure; set(gcf,'Position',[0 300 1700 1000]);
for segments = 1:17
    %define marker fill by segment
    %basal segments
    if segments == 1||segments == 2||segments == 3||segments == 4||segments == 5||segments == 6
        markeredge = "k";
        markerfill = "w";
        %mid segments
    elseif segments == 7||segments == 8||segments == 9||segments == 10||segments == 11||segments == 12
        markeredge = [0.5 0.5 0.5];
        markerfill = [0.5 0.5 0.5];
        %apical segments
    else
        markeredge = "k";
        markerfill = "k";
    end

    %figure; set(gcf,'Position',[0 300 1700 1000]);
     hold all; %set(gca,'colororder',turbo(8)); hold all;
    plot([-0.2 0.02],[-0.2 0.02],'--','Color',[0.9 0.9 0.9])
    plot(peakstrain(segments,pats),peakRSCT(segments,pats),...
        markershape,'MarkerFaceColor',markerfill,'MarkerEdgeColor',markeredge)
    if segments == 17
        text(0,-0.14,['R = ',num2str(r_rsct_pats(pats),'%.2f')],'HorizontalAlignment','right')
        if p_rsct_pats(pats) > 0.01
            text(0,-0.16,['p = ',num2str(p_rsct_pats(pats),'%.2f')],'HorizontalAlignment','right')
        else
            text(0,-0.16,'p < 0.01','HorizontalAlignment','right')
        end
        text(0,-0.18,['y = ',num2str(rsct_fit_pats(1,pats),'%.2f'),'x + ',num2str(rsct_fit_pats(2,pats),'%.2f')],'HorizontalAlignment','right')
    end
    axis([-0.2 0.02 -0.2 0.02],'square')
    xlabel('PSM Segment Strain ()')
    ylabel('CT Segment Strain (unitless)')
    %title('P_L_H_CSA')

    
end

if numel(pats) > 1
    legend('','Pat 1','Pat 2','Pat 3','Pat 4','Pat 5','Pat 6','Pat 7','Pat 8','position',[0.95 0.69 0 0])
end


%sgtitle('Agreement between Segmental PSM Work and Segmental LV Work Estimates, Patient 5')

set(findall(gcf,'-property','FontSize'),'FontSize',22)
set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',10)
set(findall(gcf,'-property','Box'),'Box','on')

%%
%%% Fisher r to z transformation %%%
% 2 = generic LVP work
% 1 = generic LVP scaled to patient-specific ES work
x1 = segMWCT_genericLVP_array;
z2 = atanh(r_MWCT_genericLVP);
z1 = atanh(r_MWCT_genericLVP_patES);
n = numel(x1);
z = (z2 - z1)./sqrt((1./(n-3))+(1./(n-3)));
%p-value
p = (1-normcdf(abs(z),0,1))*2;

error2 = 1/(sqrt(numel(x1))-3);
error1 = 1/(sqrt(numel(x1))-3);

%%% confidence interval test %%%
alpha = 0.05;
ci_l = z - 1.96*error1;
ci_u = z + 1.96*error1;

r_l = tanh(ci_l);
r_u = tanh(ci_u);

