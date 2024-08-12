%% 
% lines 3 - 73 contribute towards calculating myocardial work 
clear; clc

% NOTE: make sure your home path is the Simplified Work Estimates directory
homepath = '/Users/amandacraine/Documents/ContijochLab/repos/PSM-Work-Mapping/Simplified-Work-Estimates';
addpath([homepath,'/data/'])
addpath([homepath,'/scripts'])
addpath([homepath,'/LV Geometric Models'])
cd(homepath)
savepath = [homepath '/data-collected/'];
load RS_CT.mat
load WorkPSMLBBBCRT.mat
load seg_strain_allpats.mat

%add normality calculation at the bottom and just run stats on each
%estimate instead of saving these results
r = zeros(8,3); 
p = zeros(8,3); 
m = zeros(8,3);
b = zeros(8,3);
normality = zeros(8,3);
%% Evaluate Agreement between "MWCT" and PSM derived work
segWork_allpats = zeros(17,8);
segMWCT_allpats = zeros(17,8);
mean_Work = zeros(1,8);
segRSCTpeak_allpats = zeros(17,8);
segStrainPeak_allpats = zeros(17,8);
segRSCT_allpats = cell(1,8);
meanMWCT = zeros(1,8);
MWCT_tris = cell(1,8);
for pat = 1:8 %pat = 1:8
    foldpath = [homepath '/LV Geometric Models/BiV1-8/BiV',num2str(pat)];
    
    %work selection: work{i} or work{i}.*detj{i}, in any direction
    workPSM = work{pat}.*detj{pat};

    %load in endocardial mesh info
    modellist = dir([foldpath,'/*.obj']);
    model = readObj([foldpath,'/',modellist(1).name]);
    triangulation = model.f.v; %faces
    areas = zeros(size(triangulation,1),length(modellist));

    %pressure selection: LVP_cath, or LVPshift2
    LVP = LVP_cath;
    
    % find the mean strain for each element 
    % to do that, we need to know which faces belong to which element
    elemslist = dir([foldpath,'/*.txt']);
    elems = readtable([elemslist(1).folder,'/',elemslist(1).name]);
    ElemList = table2array(elems(:,1));

    [segMWCT,segRSCT,~,meanwork] = calculateSegmentalWorkEstimates(RS_CT{pat},LVP{pat},[],[],elementIndices,ElemList,triangulation,17,t{pat});
    segWork = calculateSegmentalWorkEstimates([],workPSM,[],[],elementIndices,ElemList,triangulation,17,t{pat});
   
    segWork_allpats(:,pat) = segWork;
    segMWCT_allpats(:,pat) = segMWCT;
    segRSCTpeak_allpats(:,pat) = min(segRSCT,[],2);
    segRSCT_allpats{pat} = segRSCT;

    mean_Work(pat) = mean(workPSM);
    meanMWCT(pat) = meanwork;

    %find peak segmental strain from PSM
    segStrainPeak_allpats(:,pat) = min(seg_strain_allpats{pat},[],2);

    for j = 1:length(triangulation)
        worktri(j) = PolyAreaSigned(RS_CT{pat}(j,:),LVP{pat});
    end
    MWCT_tris{pat} = worktri;

    
end

 save([savepath,'all_seg_work_all_pats.mat'],'segWork_allpats','segMWCT_allpats')
 save([savepath, 'all_seg_strain_all_pats.mat'],'segRSCT_allpats','seg_strain_allpats')

save([savepath 'MWCT_tris.mat'],'MWCT_tris')

 %% Plot all segmental strain and work comparisons

 figure;hold all
 plot(segWork_allpats,segMWCT_allpats,'b.');
 plot([-1 3], [-1 3],'--','Color',[0.9 0.9 0.9])
 psm = reshape(segWork_allpats,numel(segWork_allpats),1);
 test = reshape(segMWCT_allpats,numel(segMWCT_allpats),1);
 r = corr(psm,test,'type','Spearman');
 axis([-1 3 -1 3],'square')
 xlabel('Total Stress-Strain Area via PSM')
 ylabel('Pressure-Strain Area via CT')

 dum=0;
for cutoff=[min(segWork_allpats,[],'all'):0.01:max(segWork_allpats,[],'all')]
    dum=dum+1;
    segWork_label=psm>cutoff;
    classNames=[1];
    rocmets = rocmetrics(segWork_label,test,classNames);

    rocmets=addMetrics(rocmets,'ACCU');

    max_accu_idx=min(find(table2array(rocmets.Metrics(:,5)==max(rocmets.Metrics(:,5)))));
    max_accu_cutoff(dum)=table2array(rocmets.Metrics(max_accu_idx,2));
    max_accu_val(dum)=table2array(rocmets.Metrics(max_accu_idx,5));
end

cutoff=[min(segWork_allpats,[],'all'):0.01:max(segWork_allpats,[],'all')];
figure; hold all
plot(segWork_allpats,segMWCT_allpats,'.r');
plot(cutoff,max_accu_cutoff,'-k')


 figure; hold all
 plot(segStrainPeak_allpats,segRSCTpeak_allpats,'b.');
 plot([-0.2 0.1], [-0.2 0.1],'--','Color',[0.9 0.9 0.9])
 psm = reshape(segStrainPeak_allpats,numel(segStrainPeak_allpats),1);
 test = reshape(segRSCTpeak_allpats,numel(segRSCTpeak_allpats),1);
 r = corr(psm,test,'type','Spearman');
 axis([-0.2 0.1 -0.2 0.1],'square')
 xlabel('Total Peak Strain via PSM')
 ylabel('Peak Strain via CT')


  dum=0;
for cutoff=[min(segStrainPeak_allpats,[],'all'):0.001:max(segStrainPeak_allpats,[],'all')]
    dum=dum+1;
    segWork_label=psm>cutoff;
    classNames=[1];
    rocmets = rocmetrics(segWork_label,test,classNames);

    rocmets=addMetrics(rocmets,'ACCU');

    max_accu_idx=min(find(table2array(rocmets.Metrics(:,5)==max(rocmets.Metrics(:,5)))));
    max_accu_cutoff(dum)=table2array(rocmets.Metrics(max_accu_idx,2));
    max_accu_val(dum)=table2array(rocmets.Metrics(max_accu_idx,5));
end

cutoff=[min(segStrainPeak_allpats,[],'all'):0.001:max(segStrainPeak_allpats,[],'all')];
figure; hold all
plot(segStrainPeak_allpats,segRSCTpeak_allpats,'.r');
plot(cutoff(1:184),max_accu_cutoff,'-k')

%%
figure; hold all
for i = 1:8
plot(mean(seg_strain_allpats{i}),'b')
plot(mean(segRSCT_allpats{i}),'r')
end

%%
for j = 1:8
    meanMWCT(j) = mean(MWCT_tris{j});
    mean_workPSM(j) = mean(work{j}.*detj{j});
end
figure; 
hold all
plot(mean_workPSM,meanMWCT,'b.'); 
plot(mean_workPSM,mean_workPSM,'--','Color',[0.9 0.9 0.9])

