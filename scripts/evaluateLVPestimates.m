% this script analyzes how different LV pressure estimates alter our
% estimation of MW

%lines 6 - 234 collect the segmental work estimates using generic pressure
%as a surrogate for stress

clear; clc
%cd('/Users/amandacraine/Documents/ContijochLab/repos/ac-biv-mwct-validation')
%estimate 1: LV cath pressure
homepath = '/Users/amandacraine/Documents/ContijochLab/repos/CRT-PLOS-Submission-pre-repo';
cd(homepath)
savepath = [homepath '/data-collected/'];

load WorkPSMLBBBCRT.mat
load MWCT_tris.mat
load RS_CT.mat
addpath('Generic Pressure-Strain Estimate/')

%%
%addpath("readObj/")
for pat = 1:8
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
    % if LVP{pat} == LVP_cath{pat}
    %     LVP{pat} = repmat(LVP{pat},length(triangulation),1);
    % end

    % find the mean strain for each element
    % to do that, we need to know which faces belong to which element
    elemslist = dir([foldpath,'/*.txt']);
    elems = readtable([elemslist(1).folder,'/',elemslist(1).name]);
    ElemList = table2array(elems(:,1));

    %get segmental ground truth work
    segWork{pat} = calculateSegmentalWorkEstimates([],workPSM,[],[],elementIndices,ElemList,triangulation,17,t{pat});
    %get segmental MWCT
    [segMWCT{pat},segRSCT,~,meanwork{pat}] = calculateSegmentalWorkEstimates(RS_CT{pat},LVP{pat},[],[],elementIndices,ElemList,triangulation,17,t{pat});


    
    %estimate 2: completely generic LVP
    %pressures
    %paramNum = 1;
    LVPData = load('data/Generic Pressure-Strain Estimate/generic_LVP.csv');
    %LVPData(end,:) = [];
    %interpolate generic waveform to be the same size as the timing for
    %each patient
    [LVPtime,unique_time_indx,~] = unique(LVPData(:,1));
    LVP = LVPData(unique_time_indx,2);
    %make timing in terms of %RR...for now
    if LVPtime(end) > 110
        LVPtime_RR = 100.*(LVPtime./LVPtime(end));
        genericLVP_time{pat} = linspace(0,LVPtime_RR(end),numel(t{pat}));
        genericLVP{pat} = interp1(LVPtime_RR,LVP,genericLVP_time{pat})./mmHg;
    else
        genericLVP_time{pat} = linspace(0,LVPtime(end),numel(t{pat}));
        genericLVP{pat} = interp1(LVPtime,LVP,genericLVP_time{pat})./mmHg;
    end
    

    % figure; plot(LVPtime,LVP)
    % 
    % figure; plot(genericLVP_time,genericLVP{pat})


    %LVDataParam1{pat} = LVData(t{pat}+1)'./mmHg;
    % parameter 1
    for j = 1:length(triangulation)
        worktri1(j) = PolyAreaSigned(RS_CT{pat}(j,:),genericLVP{pat});
    end
    MWCT_tris1{pat} = worktri1;

    [segMWCT_genericLVP{pat},~,~,meanwork1{pat}] = calculateSegmentalWorkEstimates(RS_CT{pat},genericLVP{pat},[],[],...
        elementIndices,ElemList,triangulation,17,t{pat});
    % 
    % %estimate 3: generic LVP waveform, personalized EDP and ESP
    % % parameter 13
    % paramNum = 13;

    %scaled to patient-specific EDP
    [ESP(pat),ESP_indx(pat)] = max(LVP_cath{pat});
    genericLVP_ESP{pat} = rescale(genericLVP{pat},min(genericLVP{pat}),ESP(pat));
    % 
    for j = 1:length(triangulation)
        worktri13(j) = PolyAreaSigned(RS_CT{pat}(j,:),genericLVP_ESP{pat});
    end
    MWCT_tris13{pat} = worktri13;

    [segMWCT_genericLVP_patES{pat},~,~,meanwork13{pat}] = calculateSegmentalWorkEstimates(RS_CT{pat},genericLVP_ESP{pat},[],[],...
        elementIndices,ElemList,triangulation,17,t{pat});

    % %estimate 4: personalized timing
    % % parameter 11
    QRS = [156, 148, 162, 130, 182, 119, 140, 120]; %ms
    %if we assume a heart beats 60bpm, then one generic beat is 1000 ms
    beat_time = 1000; %ms

    %find the %of the RR interval that would be part of the QRS
    QRS_fraction = 100*(QRS(pat)./beat_time);

    %start the QRS at max dP/dt
    %can this be simplified to say that it starts at the beginning
    dPdt{pat} = diff(genericLVP{pat})./diff(genericLVP_time{pat});
    [maxdPdt,maxdPdt_indx] = max(dPdt{pat});

    %this has something to do with the correlation
    wideQRS_start = maxdPdt_indx + 1;

    %peak pressure
    [peakLVP,peakLVP_indx] = max(genericLVP{pat});

    %then stretch the curve s.t. the end of QRS ends when EDP occurs
    RS_fraction = QRS_fraction./2;
    wideQRSend_time = genericLVP_time{pat}(wideQRS_start) + RS_fraction;
    dif = abs(genericLVP_time{pat} - wideQRSend_time);
    wideQRS_end = find(dif == min(dif));

    %when does aortic valve open? we only want to know during increased LVP
    %set our curve to achieve generic LVP = 80mmHg at the end of the QRS
    genericAVO = 80./mmHg;
    %find where the pressure curve first acheives 80mmHg
    dif2 = abs(genericLVP{pat}(1:peakLVP_indx) - genericAVO);
    EDP_indx = find(dif2 == min(dif2));

    tdiff = wideQRS_end - EDP_indx;
    %genericLVP_wide = genericLVP{pat}
    alpha = wideQRS_end./EDP_indx;
    new_time = genericLVP_time{pat}.*alpha;
    %new_time = genericLVP_time{pat} + genericLVP_time{pat}(tdiff);

    %cut off unwanted x and y (x>100)
    %interp(truncated x, truncated y, linspace(0,100,215))


    %note 1/26/24: figure out what to do if the qrs adjusted time shrinks
    %the timing (time < 100)
    if new_time(end) > 100
        new_timeRR = min(find(new_time > 100));
        %wideQRS_time{pat} = new_time(1:new_timeRR);
        wideQRS_LVP{pat} = interp1(new_time(1:new_timeRR),genericLVP{pat}(1:new_timeRR),linspace(0,100,numel(t{pat})));
    else
        % %find the avg sampling rate of the last five points
        % avg_rate_change_time = mean(diff(new_time(end-5:end))); %./diff(new_time(end-4:end)))
        % %apply this sampling rate to the end of new_time until you reach
        % %100
        % %instead of a constant number, try and recirculate the LVP
        % num_added_samples = round((101 - new_time(end))./avg_rate_change_time);
        % for n = 1:num_added_samples
        %     narray_time(n) = new_time(end) + n.*avg_rate_change_time;
        %     narray_pressure(n) = genericLVP{pat}(end);
        % end
        % new_time_RR = [new_time narray_time];
        % new_pressure_RR = [genericLVP{pat} narray_pressure];
        % wideQRS_LVP{pat} = interp1(new_time_RR,new_pressure_RR,linspace(0,100,numel(t{pat})));

        %  repeat the pressure curve
        avg_rate_change_time = mean(diff(new_time(end-5:end))); %./diff(new_time(end-4:end)))
        % %apply this sampling rate to the end of new_time until you reach
        % %100
        num_added_samples = round((101 - new_time(end))./avg_rate_change_time);
        for n = 1:num_added_samples
            narray_time(n) = new_time(end) + n.*avg_rate_change_time;
        end
        new_time_RR = [new_time narray_time];
        new_pressure_RR = [genericLVP{pat} genericLVP{pat}(1:num_added_samples)];

        wideQRS_time{pat} = new_time_RR;
        wideQRS_LVP{pat} = interp1(new_time_RR,new_pressure_RR,linspace(0,100,numel(t{pat})));
        clear narray_time
    end

    

    % figure; plot(genericLVP_time{pat},genericLVP{pat});
    % hold on; plot(linspace(0,100,numel(t{pat})),wideQRS_LVP{pat})

    % %set the end of QRS to se the value of min dP/dt
    % [mindPdt,mindPdt_indx] = min(dPdt{pat});
    % %find the timing index of the end of the expanded QRS
    % wideQRSend_val = genericLVP_time{pat}(wideQRS_start) + QRS_fraction(pat);
    % dif = abs(genericLVP_time{pat} - wideQRSend_val);
    % wideQRS_end = find(dif == min(dif));

    %paramNum = 11;
    % load([LVP_paths,num2str(pat),'/Parameter',num2str(paramNum),'/LVData.mat']);
    % LVDataParam11{pat} = LVPData(t{pat}+1)'./mmHg;

    for j = 1:length(triangulation)
        worktri11(j) = PolyAreaSigned(RS_CT{pat}(j,:),wideQRS_LVP{pat});
    end
    MWCT_tris11{pat} = worktri11;

    [segMWCT11{pat},~,~,meanwork11{pat}] = calculateSegmentalWorkEstimates(RS_CT{pat},wideQRS_LVP{pat},[],[],...
        elementIndices,ElemList,triangulation,17,t{pat});

    % estimate 5: boxy pressure
    %find timing when LVP achieves ESP/2
    dif3 = abs(genericLVP{pat} - peakLVP./2);
    %dif4 = abs(genericLVP{pat}(peakLVP_indx:end) - peakLVP./2);
    halfESP1_indx = find(dif3 == min(dif3));
    halfESP2_indx = find(dif3([1:halfESP1_indx-1,halfESP1_indx+1:end]) == min(dif3([1:halfESP1_indx-1,halfESP1_indx+1:end])));

    LVPbox{pat}(1:halfESP1_indx-1) = min(genericLVP{pat});
    LVPbox{pat}(halfESP1_indx:halfESP2_indx) = peakLVP;
    LVPbox{pat}(halfESP2_indx+1:numel(t{pat})) = min(genericLVP{pat});

    % figure; hold all
    % plot(genericLVP_time{pat},genericLVP{pat})
    % plot(genericLVP_time{pat},LVPbox{pat})
    % yline(peakLVP,'k--')
    % yline(peakLVP./2,'k--')

    for j = 1:length(triangulation)
        worktri_box(j) = PolyAreaSigned(RS_CT{pat}(j,:),LVPbox{pat});
    end
    MWCT_tris_box{pat} = worktri_box;

    [segMWCT_box{pat},~,~,meanwork_box{pat}] = calculateSegmentalWorkEstimates(RS_CT{pat},LVPbox{pat},[],[],...
        elementIndices,ElemList,triangulation,17,t{pat});
    

end

segMWCT_genericLVP = cell2mat(segMWCT_genericLVP);
segMWCT_genericLVP_patES = cell2mat(segMWCT_genericLVP_patES);

save([savepath 'segMWCT_genericLVP_allpats.mat'],'segMWCT_genericLVP','segMWCT_genericLVP_patES')
save([savepath 'genericLVPs.mat'],'genericLVP','genericLVP_ESP')
%% Comparing pressure waveforms
t_indx = 1:numel(t{pat});
figure; set(gcf,'Position',[500 300 1000 700]); hold all;
pat = 1;
%cath pressure
plot(t_indx./numel(t{pat}).*100,LVP_cath{pat});
%param 1
plot(genericLVP_time{pat},genericLVP{pat});
%param 13
plot(genericLVP_time{pat},genericLVP_ESP{pat});
% %param11
% plot(linspace(0,100,numel(t{pat})),wideQRS_LVP{pat});
% %plot(wideQRS_time{pat},wideQRS_LVP{pat});
% %box estimate
% plot(genericLVP_time{pat},LVPbox{pat});


% %80 mmHg
% yline(80./mmHg,'--','Color', [0.8 0.8 0.8])
% xline(genericLVP_time{pat}(wideQRS_end))

xlabel('% RR')
ylabel('LV Pressure (kPa)')
legend('Cath LVP','Generic LVP','Patient-Specific peak LVP','Patient-Specific QRS length'),%'Generic LBBB','Generic Timing')
title(['LV Pressures for Patient ',num2str(pat)])
set(findall(gcf,'-property','FontSize'),'FontSize',28)
set(findall(gcf,'-property','LineWidth'),'LineWidth',4)
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',45)


%% Comparing Segmental MWCT estimates
segMWCTArray = cell2mat(segMWCT);
segMWCTVector = reshape(segMWCTArray,numel(segMWCTArray),1);
segMWCTArray1 = cell2mat(segMWCT_genericLVP);
segMWCTVector1 = reshape(segMWCTArray1,numel(segMWCTArray1),1);
segMWCTArray13 = cell2mat(segMWCT_genericLVP_patES);
segMWCTVector13 = reshape(segMWCTArray13,numel(segMWCTArray13),1);
segMWCTArray11 = cell2mat(segMWCT11);
segMWCTVector11 = reshape(segMWCTArray11,numel(segMWCTArray11),1);
segMWCTArray_box = cell2mat(segMWCT_box);
segMWCTVector_box = reshape(segMWCTArray_box,numel(segMWCTArray_box),1);
ytext = 2.75;
xtext = -0.75;
ytext2 = 2.45;

figure; set(gcf,'Position',[200 100 1600 500])
for pat = 1:8
    hold all;
   subplot(1,4,1); hold all;
    set(gca,'colororder',turbo(8));
    plot(segMWCT{pat},segMWCT_genericLVP{pat},'.')
   
    % norm = normalitytest(segMWCT1{pat}');
    % if norm(7,3) == 1
    %     r1 = corr(segMWCT{pat},segMWCT1{pat});
    % else
    %     r1 = corr(segMWCT{pat},segMWCT1{pat},'type','Spearman');
    % end
    if pat == 8
        plot([-0.8 2.8],[-0.8 2.8],'--','Color',[0.8 0.8 0.8])
        [r1,p1] = corr(segMWCTVector,segMWCTVector1);
        text(xtext,ytext,['r^2 = ',num2str(r1.^2,'%.2f')])
        if p1 > 0.01
            text(xtext,ytext2,['p = ',num2str(p1.^2,'%.2f')])
        else
            text(xtext,ytext2,'p < 0.01');
        end
    end
    xlabel('Cath LVP-Strain Area')
    ylabel([{'Generic LVP-Strain Area'}])%;{'(Parameter 1)'}])
    %title(['Patient ',num2str(pat),' MWCT Comparison'])
    axis([-1 3 -1 3],'square')
    box on

    subplot(1,4,2); hold all;
    set(gca,'colororder',turbo(8));
    plot(segMWCT{pat},segMWCT_genericLVP_patES{pat},'.')
    % norm = normalitytest(segMWCT13{pat}');
    % if norm(7,3) == 1
    %     r13 = corr(segMWCT{pat},segMWCT13{pat});
    % else
    %     r13 = corr(segMWCT{pat},segMWCT13{pat},'type','Spearman');
    % end
    if pat == 8
        plot([-0.8 2.8],[-0.8 2.8],'--','Color',[0.8 0.8 0.8])
        [r13,p13] = corr(segMWCTVector,segMWCTVector13);
        text(xtext,ytext,['r^2 = ',num2str(r13.^2,'%.2f')])
        if p13 > 0.01
            text(xtext,ytext2,['p = ',num2str(p13.^2,'%.2f')])
        else
            text(xtext,ytext2,'p < 0.01');
        end
    end
    xlabel('Cath LVP-Strain Area')
    ylabel([{'PS ESP LVP-Strain Area'}])%;{'(Parameter 13)'}])
    %title(['Patient ',num2str(pat),' MWCT Comparison'])
    axis([-1 3 -1 3],'square')
    box on

    subplot(1,4,3); hold all;
    set(gca,'colororder',turbo(8));
    plot(segMWCT{pat},segMWCT11{pat},'.')
    % norm = normalitytest(segMWCT11{pat}');
    % if norm(7,3) == 1
    %     r11 = corr(segMWCT{pat},segMWCT11{pat});
    % else
    %     r11 = corr(segMWCT{pat},segMWCT11{pat},'type','Spearman');
    % end
    if pat == 8
        plot([-0.8 2.8],[-0.8 2.8],'--','Color',[0.8 0.8 0.8])
        [r11,p11] = corr(segMWCTVector,segMWCTVector11);
        text(xtext,ytext,['r^2 = ',num2str(r11.^2,'%.2f')])
        if p11 > 0.01
            text(xtext,ytext2,['p = ',num2str(p11.^2,'%.2f')])
        else
            text(xtext,ytext2,'p < 0.01');
        end
    end
    xlabel('Cath LVP-Strain Area')
    ylabel({'PS QRS length LVP-Strain Area'});%{'(Parameter 11)'}])
    %title(['Patient ',num2str(pat),' MWCT Comparison'])
    axis([-1 3 -1 3],'square')
    box on

    subplot(1,4,4); hold all;
    set(gca,'colororder',turbo(8));
    plot(segMWCT{pat},segMWCT_box{pat},'.')
    % norm = normalitytest(segMWCT11{pat}');
    % if norm(7,3) == 1
    %     r11 = corr(segMWCT{pat},segMWCT11{pat});
    % else
    %     r11 = corr(segMWCT{pat},segMWCT11{pat},'type','Spearman');
    % end
    if pat == 8
        plot([-0.8 2.8],[-0.8 2.8],'--','Color',[0.8 0.8 0.8])
        [r_box,p_box] = corr(segMWCTVector,segMWCTVector_box);
        text(xtext,ytext,['r^2 = ',num2str(r_box.^2,'%.2f')])
        if p_box > 0.01
            text(xtext,ytext2,['p = ',num2str(p_box.^2,'%.2f')])
        else
            text(xtext,ytext2,'p < 0.01');
        end
    end
    xlabel('Cath LVP-Strain Area')
    ylabel({'Box LVP-Strain Area'});%{'(Parameter 11)'}])
    %title(['Patient ',num2str(pat),' MWCT Comparison'])
    axis([-1 3 -1 3],'square')
    box on

end
sgtitle('Agreement between Segmental Work Estimates with Cath LV Pressure and Generic LV Pressures')
legend('Patient 1', 'Patient 2', 'Patient 3', 'Patient 4', 'Patient 5', 'Patient 6', ...
    'Patient 7', 'Patient 8','position',[0.956 0.55 0 0])
set(findall(gcf,'-property','FontSize'),'FontSize',24)
set(findall(gcf,'-property','LineWidth'),'LineWidth',4)
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',38)
%% Comparing Segmental MWCT estimates against ground truth
segWorkArray = cell2mat(segWork);
segWorkVector = reshape(segWorkArray,numel(segWorkArray),1);
segMWCTArray = cell2mat(segMWCT);
segMWCTVector = reshape(segMWCTArray,numel(segMWCTArray),1);
segMWCTArray1 = cell2mat(segMWCT_genericLVP);
segMWCTVector1 = reshape(segMWCTArray1,numel(segMWCTArray1),1);
segMWCTArray13 = cell2mat(segMWCT_genericLVP_patES);
segMWCTVector13 = reshape(segMWCTArray13,numel(segMWCTArray13),1);
segMWCTArray11 = cell2mat(segMWCT11);
segMWCTVector11 = reshape(segMWCTArray11,numel(segMWCTArray11),1);
segMWCTArray_box = cell2mat(segMWCT_box);
segMWCTVector_box = reshape(segMWCTArray_box,numel(segMWCTArray_box),1);

ytext = 2.75;
ytext2 = 2.35;
xtext = -0.75; 
figure; set(gcf,'Position',[200 200 1600 400])
for pat = 1:8
    subplot(1,3,1); hold all
    set(gca,'colororder',turbo(8));
    plot(segWork{pat},segMWCT{pat},'.')
    % norm = normalitytest(segMWCT{pat}');
    % if norm(7,3) == 1
    %     r = corr(segWork{pat},segMWCT{pat});
    % else
    %     r = corr(segWork{pat},segMWCT{pat},'type','Spearman');
    % end
    if pat == 8
        plot([-0.8 2.8],[-0.8 2.8],'--','Color',[0.8 0.8 0.8])
        [r,p] = corr(segWorkVector,segMWCTVector);
        %text(xtext,ytext,['r^2 = ',num2str(r.^2,'%.2f')])
        text(xtext,ytext,['r = ',num2str(r,'%.2f')])
        if p > 0.01
            text(xtext,ytext2,['p = ',num2str(p,'%.2f')])
        else
            text(xtext,ytext2,'p < 0.01')
        end
    end
    xlabel('Ground Truth')
    ylabel('Cath LVP-Strain Area')
    %title(['Patient ',num2str(pat),' MWCT Comparison'])
    axis([-1 3 -1 3],'square')
    box on

    subplot(1,3,2); hold all;
    set(gca,'colororder',turbo(8));
    plot(segWork{pat},segMWCT_genericLVP{pat},'.')
    % norm = normalitytest(segMWCT1{pat}');
    % if norm(7,3) == 1
    %     r1 = corr(segWork{pat},segMWCT1{pat});
    % else
    %     r1 = corr(segWork{pat},segMWCT1{pat},'type','Spearman');
    % end
    if pat == 8
        plot([-0.8 2.8],[-0.8 2.8],'--','Color',[0.8 0.8 0.8])
        [r1,p1] = corr(segWorkVector,segMWCTVector1);
        %text(xtext,ytext,['r^2 = ',num2str(r1.^2,'%.2f')])
        text(xtext,ytext,['r = ',num2str(r1,'%.2f')])
        if p1 > 0.01
            text(xtext,ytext2,['p = ',num2str(p1,'%.2f')])
        else
            text(xtext,ytext2,'p < 0.01')
        end
    end
    xlabel('Ground Truth')
    ylabel([{'Generic LVP-Strain Area'}]);%;{'(Parameter 1)'}])
    %title(['Patient ',num2str(pat),' MWCT Comparison'])
    axis([-1 3 -1 3],'square')
    box on

    subplot(1,3,3); hold all;
    set(gca,'colororder',turbo(8));
    plot(segWork{pat},segMWCT_genericLVP_patES{pat},'.')
    % norm = normalitytest(segMWCT13{pat}');
    % if norm(7,3) == 1
    %     r13 = corr(segWork{pat},segMWCT13{pat});
    % else
    %     r13 = corr(segWork{pat},segMWCT13{pat},'type','Spearman');
    % end
    if pat == 8
        plot([-0.8 2.8],[-0.8 2.8],'--','Color',[0.8 0.8 0.8])
        [r13,p13] = corr(segWorkVector,segMWCTVector13);
        %text(xtext,ytext,['r^2 = ',num2str(r13.^2,'%.2f')])
        text(xtext,ytext,['r = ',num2str(r13,'%.2f')])
        if p13 > 0.01
            text(xtext,ytext2,['p = ',num2str(p13,'%.2f')])
        else
            text(xtext,ytext2,'p < 0.01')
        end
    end
    xlabel('Ground Truth')
    ylabel([{'Patient-Specific '};{'peak LVP-Strain Area'}]);%{'(Parameter 13)'}])
    %title(['Patient ',num2str(pat),' MWCT Comparison'])
    axis([-1 3 -1 3],'square')
    box on
    % 
    % subplot(1,5,4); hold all
    % set(gca,'colororder',turbo(8));
    % plot(segWork{pat},segMWCT11{pat},'.')
    % % norm = normalitytest(segMWCT11{pat}');
    % % if norm(7,3) == 1
    % %     r11 = corr(segWork{pat},segMWCT11{pat});
    % % else
    % %     r11 = corr(segWork{pat},segMWCT11{pat},'type','Spearman');
    % % end
    % if pat == 8
    %     plot([-0.8 2.8],[-0.8 2.8],'--','Color',[0.8 0.8 0.8])
    %     [r11,p11] = corr(segWorkVector,segMWCTVector11);
    %     text(xtext,ytext,['r^2 = ',num2str(r11.^2,'%.2f')])
    %     if p11 > 0.01
    %         text(xtext,ytext2,['p = ',num2str(p11,'%.2f')])
    %     else
    %         text(xtext,ytext2,'p < 0.01')
    %     end
    % end
    % xlabel('Ground Truth')
    % ylabel([{'PS QRS Timing Timing LVP-Strain Area'}]);%{'(Parameter 11)'}])
    % %title(['Patient ',num2str(pat),' MWCT Comparison'])
    % axis([-1 3 -1 3],'square')
    % box on
    % 
    % subplot(1,5,5); hold all
    % set(gca,'colororder',turbo(8));
    % plot(segWork{pat},segMWCT_box{pat},'.')
    % % norm = normalitytest(segMWCT11{pat}');
    % % if norm(7,3) == 1
    % %     r11 = corr(segWork{pat},segMWCT11{pat});
    % % else
    % %     r11 = corr(segWork{pat},segMWCT11{pat},'type','Spearman');
    % % end
    % if pat == 8
    %     plot([-0.8 2.8],[-0.8 2.8],'--','Color',[0.8 0.8 0.8])
    %     [r_box,p_box] = corr(segWorkVector,segMWCTVector_box);
    %     text(xtext,ytext,['r^2 = ',num2str(r_box.^2,'%.2f')])
    %     if p11 > 0.01
    %         text(xtext,ytext2,['p = ',num2str(p_box,'%.2f')])
    %     else
    %         text(xtext,ytext2,'p < 0.01')
    %     end
    % end
    % xlabel('Ground Truth')
    % ylabel([{'Box LVP-Strain Area'}]);%{'(Parameter 11)'}])
    % %title(['Patient ',num2str(pat),' MWCT Comparison'])
    % axis([-1 3 -1 3],'square')
    % box on

end
sgtitle('Agreement between Segmental Work Estimates and Ground Truth Segmental Work')
legend('Patient 1', 'Patient 2', 'Patient 3', 'Patient 4', 'Patient 5', 'Patient 6', ...
    'Patient 7', 'Patient 8','position',[0.956 0.51 0 0])
set(findall(gcf,'-property','FontSize'),'FontSize',22)
set(findall(gcf,'-property','LineWidth'),'LineWidth',4)
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',38)
