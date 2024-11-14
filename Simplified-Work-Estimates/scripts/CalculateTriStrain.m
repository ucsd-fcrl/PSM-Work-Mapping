%Lines 4 - 26 are calculating RSCT.
% lines 27 - 170 collects seg_strain_allpats

clear; clc
% NOTE: make sure your home path is the Simplified Work Estimates directory
homepath = '/Users/amandacraine/Documents/ContijochLab/repos/PSM-Work-Mapping/Simplified-Work-Estimates';
addpath([homepath,'/data/'])
addpath([homepath,'/scripts'])
addpath([homepath,'/LV Geometric Models'])
cd(homepath)
savepath = [homepath '/data-collected/'];
%load in the sample files
RS_CT = cell(1,8);
for pat = 1:8
    foldpath = [homepath '/LV Geometric Models/BiV1-8/BiV',num2str(pat)];    
    modellist = dir([foldpath,'/*.obj']);
    %addpath('readObj/')
    
    model = readObj([foldpath,'/',modellist(1).name]);
    triangulation = model.f.v; %faces
    areas = zeros(size(triangulation,1),length(modellist));
    
    RS_CT{pat} = calculate_RSCT(foldpath,modellist,triangulation,areas);

end
save([savepath 'RS_CT.mat'],'RS_CT')    
 %% Big Scatter Plot comparing strains
%  clear; clc
% addpath('readObj/')

load RS_CT.mat
load dataPSMLBBBCRT.mat

r = zeros(8,4);
p = zeros(8,4);
m = zeros(8,4);
b = zeros(8,4);

ES_time = [9140,9200,9168,9108,10580,9204,9072,9132];

seg_strain_allpats = cell(1,8);
for pat = 1:8
    
    foldpath = ['/Users/amandacraine/Documents/ContijochLab/BiVFunctionModels/BiV1-8/BiV',num2str(pat)];

    modellist = dir([foldpath,'/*.obj']);
    model = readObj([foldpath,'/',modellist(1).name]);
    triangulation = model.f.v; %faces
    areas = zeros(size(triangulation,1),length(modellist));

    % find the mean strain for each element
    % to do that, we need to know which faces belong to which element
    elemslist = dir([foldpath,'/*.txt']);
    elems = readtable([elemslist(1).folder,'/',elemslist(1).name]);
    ElemList = table2array(elems(:,1));

    elemnames = cell(length(elemslist),1);
    for k = 1:length(elemslist)
        elemnames{k} = elemslist(k).name;
        %ESindx = find(strcmp(elemnames,'LVPoints913200.txt'));
    end
    ESindx = find(strcmp(elemnames,['LVPoints',num2str(ES_time(pat)),'00.txt']));
   

    Elem2Tri = zeros(size(triangulation));
    for i = 1:size(triangulation,1)
        for j = 1:size(triangulation,2)
            Elem2Tri(i,j) = ElemList(triangulation(i,j));
        end
    end

    segmentElems = cell(17,1);
    segmentElems{2} = [1,2];
    segmentElems{1} = [5];
    segmentElems{6} = [6];
    segmentElems{5} = [9];
    segmentElems{4} = [10];
    segmentElems{3} = [13,14];

    segmentElems{8} = [18,20];
    segmentElems{7} = [26];
    segmentElems{12} = [28];
    segmentElems{11} = [34];
    segmentElems{10} = [36];
    segmentElems{9} = [42,44];

    segmentElems{14} = [17,19,41,43];
    segmentElems{15} = [35];
    segmentElems{16} = [27,33];
    segmentElems{13} = [25];
    segmentElems{17} = [49,50,53,54,57,58,61,62];

    % if all the vertices in Elem2tri have the same elements (we assume they all do but just to check)
    % then store the face as that one element
    TriElems = zeros(length(Elem2Tri),1);
    for i = 1:length(Elem2Tri)
        if isequal(Elem2Tri(i,1),Elem2Tri(i,2),Elem2Tri(i,3))
            TriElems(i) = Elem2Tri(i,1);
        end
    end

    TriSegElems = zeros(length(triangulation),1);
    for i = 1:length(triangulation)
        for k = 1:length(segmentElems)
            if any((segmentElems{k} == TriElems(i)))
                TriSegElems(i) = k;
            end
        end
    end

    pointsSeg = zeros(length(elementIndices),1);
    for i = 1:length(elementIndices)
        for k = 1:length(segmentElems)
            if any((segmentElems{k} == elementIndices(i)))
                pointsSeg(i) = k;
            end
        end
    end
    i = pat;
    numTimePoints = length(t{i});
    segRSCT = zeros(17,numTimePoints);
    seg_strain = zeros(17,numTimePoints);
    % seg_strainpf = zeros(17,numTimePoints);
    % seg_strainpc = zeros(17,numTimePoints);
    % seg_strainps = zeros(17,numTimePoints);
    % seg_strainF = zeros(17,numTimePoints);
    % seg_strainC = zeros(17,numTimePoints);
    % seg_strainS = zeros(17,numTimePoints);

     for region = 1:17
         segRSCT(region,:) = mean(RS_CT{i}(TriSegElems == region,:));
         seg_strain(region,:) = mean(totalstrain{i}(pointsSeg == region,:));
    %     seg_strainpf(region,:) = mean(strainpf{i}(pointsSeg == region,:));
    %     seg_strainpc(region,:) = mean(strainpc{i}(pointsSeg == region,:));
    %     seg_strainps(region,:) = mean(strainps{i}(pointsSeg == region,:));
     end

    % for time = 1:size(segRSCT,2)
    %     seg_strainpf(:,time) = (seg_strainF(:,time)  - seg_strainF(:,1));%./segStrainF(:,1);
    %     seg_strainpc(:,time) = (seg_strainC(:,time)  - seg_strainC(:,1));%./segStrainF(:,1);
    %     seg_strainps(:,time) = (seg_strainS(:,time)  - seg_strainS(:,1));%./segStrainF(:,1);
    % end

    peakstrainPSM = [min(seg_strain,[],2)];%, min(seg_strainpf,[],2),min(seg_strainpc,[],2),min(seg_strainps,[],2)];
    peakstrainRSCT = min(segRSCT,[],2);
    ESstrainPSM = [seg_strain(:,ESindx)];%, seg_strainpf(:,ESindx),seg_strainpc(:,ESindx),seg_strainps(:,ESindx)];
    ES_RSCT = segRSCT(:,ESindx);
    % %straintype = {peakstrainPSM,peakstrainRSCT;ESstrainPSM,ES_RSCT};
    PSMstraintype = {peakstrainPSM;ESstrainPSM};
    RSCTstraintype = {peakstrainRSCT;ES_RSCT};
    straintype_indx = 1;
    
        estimate1 = data_analysis(PSMstraintype{straintype_indx}(:,1),RSCTstraintype{straintype_indx},1,1);%data_analysis(min(seg_strain,[],2),min(segRSCT,[],2),1);
        % estimate2 = data_analysis(PSMstraintype{straintype_indx}(:,2),RSCTstraintype{straintype_indx},1,1);
        % estimate3 = data_analysis(PSMstraintype{straintype_indx}(:,3),RSCTstraintype{straintype_indx},1,1);
        % estimate4 = data_analysis(PSMstraintype{straintype_indx}(:,4),RSCTstraintype{straintype_indx},1,1);
        % 
        normality(pat,:) = [estimate1(1)];%,estimate2(1),estimate3(1)];
        r(i,:) = [estimate1(8)]; %, estimate2(8), estimate3(8), estimate4(8)];
        p(i,:) = [estimate1(9)]; %, estimate2(9), estimate3(9),estimate4(9)];
        m(i,:) = [estimate1(10)]; %, estimate2(10), estimate3(10),estimate4(10)];
        b(i,:) = [estimate1(11)]; %, estimate2(11), estimate3(11),estimate4(11)];
        % title =  {'Peak Strain';'End-Systolic Strain'};
  % [fig] = makeDirectionalScatterPlots(PSMstraintype{straintype_indx},RSCTstraintype{straintype_indx},i,r,p,b,m,title{straintype_indx});
 % saveas(gcf,['/Users/amandacraine/Documents/ContijochLab/repos/ac-biv-mwct-validation/Results/RSCT vs Ground Truth/Pat ',num2str(i),'_',title{straintype_indx},'_comp.jpg'])

  seg_strain_allpats{pat} = seg_strain; 
end

save([savepath 'seg_strain_allpats.mat'],'seg_strain_allpats')

 r_median = median(r);
 q1_r1 = prctile(r(:,1),25);
% q1_r2 = prctile(r(:,2),25);
% q1_r3 = prctile(r(:,3),25);
% q1_r4 = prctile(r(:,4),25);
 q4_r1 = prctile(r(:,1),75);
% q4_r2 = prctile(r(:,2),75);
% q4_r3 = prctile(r(:,3),75);
% q4_r4 = prctile(r(:,4),75);
% iqr_r = [q1_r1,q1_r2,q1_r3,q1_r4;q4_r1,q4_r2,q4_r3,q4_r4];
 m_median = median(m);
 q1_m1 = prctile(m(:,1),25);
% q1_m2 = prctile(m(:,2),25);
% q1_m3 = prctile(m(:,3),25);
% q1_m4 = prctile(m(:,4),25);
 q4_m1 = prctile(m(:,1),75);
% q4_m2 = prctile(m(:,2),75);
% q4_m3 = prctile(m(:,3),75);
% q4_m4 = prctile(m(:,4),75);
% iqr_m = [q1_m1,q1_m2,q1_m3,q1_m4;q4_m1,q4_m2,q4_m3,q4_m4];

%% Plot ED and ES regional strain on the mesh
clear title
pat = 5;

foldpath = ['/Users/amandacraine/Documents/ContijochLab/BiVFunctionModels/BiV1-8/BiV',num2str(pat)];
modellist = dir([foldpath,'/*.obj']);
model = readObj([foldpath,'/',modellist(1).name]); %End-diastole
framepts = model.v; %individuals vertex pts
triangulation = model.f.v; %faces

elemslist = dir([foldpath,'/*.txt']);
elems = readtable([elemslist(1).folder,'/',elemslist(1).name]);
ElemList = table2array(elems(:,1));

elemnames = cell(length(elemslist),1);
for k = 1:length(elemslist)
    elemnames{k} = elemslist(k).name;
    %ESindx = find(strcmp(elemnames,'LVPoints913200.txt'));
end

ES_time = [9140,9200,9168,9108,10580,9204,9072,9132];
%ESindx = find(strcmp(elemnames,['LVPoints',num2str(ES_time(pat)),'00.txt']));
ESindx = find(strcmp(elemnames,['LVPoints',num2str(ES_time(pat)),'00.txt']));

wall_angle = [90 + -1.797075218781851e+02,180 + -1.797075218781851e+02, 270 + -1.797075218781851e+02, -1.797075218781851e+02];
wall_name = [{'Anterior'}, {'Lateral'},{'Posterior'},{'Septal'}];
figure; set(gcf,'Position',[200 200 1500 500])
for i = 1:4
subplot(2,4,i)
a = trisurf(triangulation,framepts(:,3),framepts(:,2),-framepts(:,1)); hold on
a.CData = RS_CT{pat}(:,1);%TriSegElems;
%c = pink; %turbo;
end_colormap1 = [145/255, 95/255, 109/255];
c1 = [linspace(1,end_colormap1(1),256)',linspace(1,end_colormap1(2),256)',linspace(1,end_colormap1(3),256)'];
end_colormap2 = [30/255 18/255 26/255];
c2 = [linspace(end_colormap1(1),end_colormap2(1),256)',linspace(end_colormap1(2),end_colormap2(2),256)',linspace(end_colormap1(3),end_colormap2(3),256)'];
c = [c1;c2];
colormap(flipud(c(1:2:end,:)))
clim([-0.1 0.02])%colorbar('Ticks',1:17,'TickLabels',["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17"])
grid off
axis off

zlim([-8 2]);
ylim([-6 4]);
xlim([-4 4]);
view([wall_angle(i),1.758832116788327])
title([wall_name{i},' Wall']);
axis square

subplot(2,4,i+4)
a = trisurf(triangulation,framepts(:,3),framepts(:,2),-framepts(:,1)); hold on
a.CData = RS_CT{pat}(:,ESindx);%TriSegElems;
%turbo;
colormap(flipud(c(1:2:end,:)))
clim([-0.1 0])%colorbar('Ticks',1:17,'TickLabels',["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17"])
grid off
axis off

zlim([-8 2]);
ylim([-6 4]);
xlim([-4 4]);
view([wall_angle(i),1.758832116788327])
title([wall_name{i},' Wall']);
axis square
end

cb = colorbar('Ticks',[-0.1:0.05:0]);
set(cb,'position',[.94 .16 0.0107 0.75])
set(findall(gcf,'-property','FontSize'),'FontSize',25)

%% Bullseye plot of RSCT

addpath(genpath('AHABullsEye'))
figure; set(gcf,'Position',[800 100 800 650])
c = createBullseye([0 0.5 1 0; 0.5 1 4 45;1 1.5 6 0; 1.5 2 6 0]);
set(c,'Color','k','LineWidth',2)
box off

for region = 1:17
    straincurve = mean(RS_CT{pat}(TriSegElems == region,:));
    %straincurve2 = mean(totalstrain{pat}(pointsSeg == region,:));
    data(region) = straincurve(ESindx);
    %data(region) = min(straincurve);
    % data1(region) = straincurve1(ESindx);
    % data2(region) = straincurve2(ESindx);

end

% Example 1 of filling the bullseye, vector by vector
fillBullseye([data(17)],0,1,0,360)
fillBullseye([data(13:16)],0.5,1,45,405);
fillBullseye([data(7:12)],1,1.5,60,420);
fillBullseye([data(1:6)],1.5,2,60,420);

%colormap(turbo)
end_colormap1 = [145/255, 95/255, 109/255];
c1 = [linspace(1,end_colormap1(1),256)',linspace(1,end_colormap1(2),256)',linspace(1,end_colormap1(3),256)'];
end_colormap2 = [30/255 18/255 26/255];
c2 = [linspace(end_colormap1(1),end_colormap2(1),256)',linspace(end_colormap1(2),end_colormap2(2),256)',linspace(end_colormap1(3),end_colormap2(3),256)'];
c = [c1;c2];
colormap(flipud(c(1:2:end,:)))
clim([-0.1 0])
%colorbar
cb = colorbar('Ticks',[-0.1:0.05:0]);
set(cb,'position',[.865 .23 0.02 0.55])

text(-0.05,1.75,num2str(data(1),'%.2f'),'HorizontalAlignment','center','Color',[0.9 0.9 0.9])
text(-1.55,0.85,num2str(data(2),'%.2f'),'HorizontalAlignment','center')
text(-1.55,-0.85,num2str(data(3),'%.2f'),'HorizontalAlignment','center')
text(-0.05,-1.75,num2str(data(4),'%.2f'),'HorizontalAlignment','center')
text(1.5,-0.85,num2str(data(5),'%.2f'),'HorizontalAlignment','center','Color',[0.9 0.9 0.9])
text(1.5,0.85,num2str(data(6),'%.2f'),'HorizontalAlignment','center','Color',[0.9 0.9 0.9])

text(-0.05,1.25,num2str(data(7),'%.2f'),'HorizontalAlignment','center')
text(-1.1,0.6,num2str(data(8),'%.2f'),'HorizontalAlignment','center')
text(-1.15,-0.6,'0','HorizontalAlignment','center')
text(0,-1.25,num2str(data(10),'%.2f'),'HorizontalAlignment','center')
text(1.1,-0.6,num2str(data(11),'%.2f'),'HorizontalAlignment','center','Color',[0.9 0.9 0.9])
text(1.1,0.6,num2str(data(12),'%.2f'),'HorizontalAlignment','center')

text(0,0.75,num2str(data(13),'%.2f'),'HorizontalAlignment','center')
text(-0.75,0,'0','HorizontalAlignment','center')
text(0,-0.75,num2str(data(15),'%.2f'),'HorizontalAlignment','center')
text(0.75,0,num2str(data(16),'%.2f'),'HorizontalAlignment','center')

text(0,0,num2str(data(17),'%.2f'),'HorizontalAlignment','center')
set(findall(gcf,'-property','FontSize'),'FontSize',20)


