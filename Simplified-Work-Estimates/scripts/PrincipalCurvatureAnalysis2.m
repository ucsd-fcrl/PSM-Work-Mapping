%This script calculates principal curvatures by fitting the geometric model
%to an ellipsoid
 
%lines 30 - 125 actually collect and store segmental MW estimates with
%segmental wall stresses and strains

% both segMWCT_effrad_allpats.mat and Laplace_measurements.mat variables
% calculated and saved here

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
load("WorkPSMLBBBCRT.mat")
savepath = [homepath '/data-collected/'];

%% test
%load in the mesh
pat = 5;
foldpath = [homepath '/LV Geometric Models/BiV1-8/BiV',num2str(pat)];
modellist = dir([foldpath,'/*.obj']);
model = readObj([foldpath,'/',modellist(1).name]); %End-diastole
framepts = model.v; %individuals vertex pts
triangulation = model.f.v; %faces

% Plot the LV mesh
figure; set(gcf,'Position',[200 200 500 500])
a = trisurf(triangulation,framepts(:,3),framepts(:,2),-framepts(:,1),'FaceColor',[0.6350 0.0780 0.1840],'FaceAlpha','0.8');
grid off
axis off
axis square

%% can we fit an arbitrary ellipsoid to a segment or patch?
clear; clc
addpath([homepath,'scripts/ellipsoid_fit/ellipsoid_fit/'])
load RS_CT.mat
load WorkPSMLBBBCRT.mat

warning off
%datapath = '/Users/amandacraine/Documents/ContijochLab/repos/ac-continuity/';
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

for pat = 1:8
    disp(['Analyzing patient ',num2str(pat)])
    foldpath = [homepath '/LV Geometric Models/BiV1-8/BiV',num2str(pat)];
    for time = 1:length(t{pat})
        modellist = dir([foldpath,'/*.obj']);
        model = readObj([foldpath,'/',modellist(time).name]); %End-diastole
        framepts = model.v; %individuals vertex pts
        triangulation = model.f.v; %faces
        thickness_data{pat} = readmatrix([homepath '/Segmental Wall Thickness Measurements/BiV',num2str(pat),'/WallThicknessWholeLVTimeVar.txt']);

        elemslist = dir([foldpath,'/*.txt']);
        elems = readtable([elemslist(1).folder,'/',elemslist(time).name]);
        ElemList = table2array(elems(:,1));

        Elem2Tri = zeros(size(triangulation));
        for i = 1:size(triangulation,1)
            for j = 1:size(triangulation,2)
                Elem2Tri(i,j) = ElemList(triangulation(i,j));
            end
        end

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

        for region = 1:17
            %disp(['Evaluating region ',num2str(region)])
            % get segmental longitudinal radius

            x = framepts(triangulation(TriSegElems == region,:),3);
            y = framepts(triangulation(TriSegElems == region,:),2);
            z = -framepts(triangulation(TriSegElems == region,:),1);

            [center, radii(:,region), evecs, v, chi2] = ellipsoid_fit([x,y,z], '0xy');

            r1 = radii(3,region);
            r2 = radii(1,region);
            r_eff(region,time) = r1.*r2./(r1+r2);

        end

    end
    r_eff_allpats{pat} = r_eff;


    segMWCTLP(:,pat) = calculateSegmentalWorkEstimates(RS_CT{pat},LVP_cath{pat},r_eff_allpats{pat}(:,1),thickness_data{pat}(:,1),elementIndices,ElemList,triangulation,17,t{pat});
    segMWCTLPTV(:,pat) = calculateSegmentalWorkEstimates(RS_CT{pat},LVP_cath{pat},r_eff_allpats{pat},thickness_data{pat},elementIndices,ElemList,triangulation,17,t{pat});
    clear r_eff
end

save([savepath 'segMWCT_effrad_allpats.mat'],'segMWCTLP','segMWCTLPTV')
save([savepath 'Laplace_measurements.mat'],'r_eff_allpats','thickness_data')

%%
load all_seg_work_all_pats.mat
load segMWCT_effrad_allpats.mat
figure; set(gcf,'Position',[200 500 1400 1000]); 

subplot(2,3,1); set(gca,'colororder',turbo(8)); hold all; 
plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
plot(segWork_allpats,segMWCT_allpats,'.');
axis([-0.5 4 -5 4],'square')
[r,p] = corr(reshape(segWork_allpats,numel(segWork_allpats),1),reshape(segMWCT_allpats,numel(segMWCT_allpats),1),'type','spearman');
fitvals = polyfit(reshape(segWork_allpats,numel(segWork_allpats),1),reshape(segMWCT_allpats,numel(segMWCT_allpats),1),1);
text(1.5,-3.2,['r = ',num2str(r,'%.2f')])
text(1.5,-4,['y = ',num2str(fitvals(1),'%.2f'),'x + ',num2str(fitvals(2),'%.2f')])
ylabel('LV Pressure-Strain Area via CT')
xlabel('Stress-Strain Area via PSM')

subplot(2,3,2); set(gca,'colororder',turbo(8)); hold all
plot([-5 12],[-5 12],'--','Color',[0.9 0.9 0.9])
plot(segWork_allpats,segMWCT_LLP_allpats,'.');
[r,p] = corr(reshape(segWork_allpats,numel(segWork_allpats),1),reshape(segMWCT_LLP_allpats,numel(segMWCT_LLP_allpats),1),'type','spearman');
fitvals = polyfit(reshape(segWork_allpats,numel(segWork_allpats),1),reshape(segMWCT_LLP_allpats,numel(segMWCT_LLP_allpats),1),1);
text(1.5,-1.6,['r = ',num2str(r,'%.2f')])
text(1.5,-3.1113,['y = ',num2str(fitvals(1),'%.2f'),'x + ',num2str(fitvals(2),'%.2f')])
axis([-0.5 4 -5 12],'square')
ylabel({'Static Laplace';'Stress-Strain Area via CT'})
xlabel('Stress-Strain Area via PSM')

subplot(2,3,3); set(gca,'colororder',turbo(8)); hold all
plot([-5 12],[-5 12],'--','Color',[0.9 0.9 0.9])
plot(segWork_allpats,segMWCT_LLP_TV_allpats,'.');
[r,p] = corr(reshape(segWork_allpats,numel(segWork_allpats),1),reshape(segMWCT_LLP_TV_allpats,numel(segMWCT_LLP_TV_allpats),1),'type','spearman');
fitvals = polyfit(reshape(segWork_allpats,numel(segWork_allpats),1),reshape(segMWCT_LLP_TV_allpats,numel(segMWCT_LLP_TV_allpats),1),1);
text(1.5,-1.6,['r = ',num2str(r,'%.2f')])
text(1.5,-3.1113,['y = ',num2str(fitvals(1),'%.2f'),'x + ',num2str(fitvals(2),'%.2f')])
axis([-0.5 4 -5 12],'square')
ylabel({'Time-Varying Laplace'; 'Stress-Strain Area via CT'})
xlabel('Stress-Strain Area via PSM')
%legend('','Pat 1','Pat 2','Pat 3','Pat 4','Pat 5','Pat 6','Pat 7','Pat 8','position',[0.95 0.5 0 0])


subplot(2,3,4); set(gca,'colororder',turbo(8)); hold all
plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
plot(segWork_allpats,segMWCT_allpats,'.');
axis([-0.5 4 -5 4],'square')
[r,p] = corr(reshape(segWork_allpats,numel(segWork_allpats),1),reshape(segMWCT_allpats,numel(segMWCT_allpats),1),'type','spearman');
fitvals = polyfit(reshape(segWork_allpats,numel(segWork_allpats),1),reshape(segMWCT_allpats,numel(segMWCT_allpats),1),1);
text(1.5,-3.2,['r = ',num2str(r,'%.2f')])
text(1.5,-4,['y = ',num2str(fitvals(1),'%.2f'),'x + ',num2str(fitvals(2),'%.2f')])
ylabel('LV Pressure-Strain Area via CT')
xlabel('Stress-Strain Area via PSM')

subplot(2,3,5); set(gca,'colororder',turbo(8)); hold all
plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
plot(segWork_allpats,segMWCTLP,'.');
[r,p] = corr(reshape(segWork_allpats,numel(segWork_allpats),1),reshape(segMWCTLP,numel(segMWCTLP),1),'type','spearman');
fitvals = polyfit(reshape(segWork_allpats,numel(segWork_allpats),1),reshape(segMWCTLP,numel(segMWCTLP),1),1);
text(1.5,-3.2,['r = ',num2str(r,'%.2f')])
text(1.5,-4,['y = ',num2str(fitvals(1),'%.2f'),'x + ',num2str(fitvals(2),'%.2f')])
axis([-0.5 4 -5 4],'square')
ylabel({'Static Laplace';'Stress-Strain Area via CT'})
xlabel('Stress-Strain Area via PSM')

subplot(2,3,6); set(gca,'colororder',turbo(8)); hold all
plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
plot(segWork_allpats,segMWCTLPTV,'.');
[r,p] = corr(reshape(segWork_allpats,numel(segWork_allpats),1),reshape(segMWCTLPTV,numel(segMWCTLPTV),1),'type','spearman');
fitvals = polyfit(reshape(segWork_allpats,numel(segWork_allpats),1),reshape(segMWCTLPTV,numel(segMWCTLPTV),1),1);
text(1.5,-3.2,['r = ',num2str(r,'%.2f')])
text(1.5,-4,['y = ',num2str(fitvals(1),'%.2f'),'x + ',num2str(fitvals(2),'%.2f')])
axis([-0.5 4 -5 4],'square')
ylabel({'Time-Varying Laplace'; 'Stress-Strain Area via CT'})
xlabel('Stress-Strain Area via PSM')
legend('','Pat 1','Pat 2','Pat 3','Pat 4','Pat 5','Pat 6','Pat 7','Pat 8','position',[0.95 0.75 0 0])

sgtitle('Segmental MW Estimates with Uniaxial Radius Measurements')
text(-4.5,5,'Segmental MW Estimates with Biaxial Radius Measurements','HorizontalAlignment','center')
set(findall(gcf,'-property','FontSize'),'FontSize',22)
set(findall(gcf,'-property','LineWidth'),'LineWidth',4)
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',38)
set(findall(gcf,'-property','Box'),'Box','on')

%%

% bland altman...for fun
diff_PSA = segWork_allpats - segMWCT_allpats;
mean_PSA = (segWork_allpats + segMWCT_allpats)./2;

diff_MWLP = segWork_allpats - segMWCTLP;
mean_MWLP = (segWork_allpats + segMWCTLP)./2;

diff_MWLPTV = segWork_allpats - segMWCTLPTV;
mean_MWLPTV = (segWork_allpats + segMWCTLPTV)./2;

figure; 
subplot(1,3,1)
plot(mean_PSA,diff_PSA,'.')
yline(median(diff_PSA,'all'),'r--')
yline(prctile(diff_PSA,2.5,'all'),'b--')
yline(prctile(diff_PSA,97.5,'all'),'b--')
xlabel('Mean MW')
ylabel('MW_P_S_M - MW_C_T')
title('Pressure-Strain Area')
axis([-3 3 -3 5])

subplot(1,3,2)
plot(mean_MWLP,diff_MWLP,'.')
yline(median(diff_MWLP,'all'),'r--')
yline(prctile(diff_MWLP,2.5,'all'),'b--')
yline(prctile(diff_MWLP,97.5,'all'),'b--')
xlabel('Mean MW')
ylabel('MW_P_S_M - MW_C_T')
title('Static Laplace Stress-Strain Area')
axis([-3 3 -3 5])

subplot(1,3,3)
plot(mean_MWLPTV,diff_MWLPTV,'.')
yline(median(diff_MWLPTV,'all'),'r--')
yline(prctile(diff_MWLPTV,2.5,'all'),'b--')
yline(prctile(diff_MWLPTV,97.5,'all'),'b--')
xlabel('Mean MW')
ylabel('MW_P_S_M - MW_C_T')
title('Time-Varying Laplace Stress-Strain Area')
axis([-3 3 -3 5])





%% Patient specific plots
load segMWCT_effrad_allpats.mat
for pat = 1
    figure; set(gcf,'Position',[200 500 1400 500])
    subplot(1,3,1); hold all
    plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
    plot(segWork_allpats(:,pat),segMWCT_allpats(:,pat),'b.');
    axis([-0.5 4 -5 4],'square')
    [r1,p] = corr(segWork_allpats(:,pat),segMWCT_allpats(:,pat),'type','spearman');
    fit1 = polyfit(segWork_allpats(:,pat),segMWCT_allpats(:,pat),1);
    text(1.5,-3.2,['r^2 = ',num2str(r1.^2,'%.2f')])
    text(1.5,-4,['y = ',num2str(fit1(1),'%.2f'),'x + ',num2str(fit1(2),'%.2f')])
    ylabel('LV Pressure-Strain Area via CT')
    xlabel('Stress-Strain Area via PSM')
    
    subplot(1,3,2); hold all
    plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
    plot(segWork_allpats(:,pat),segMWCTLP(:,pat),'b.');
    [r2,p] = corr(segWork_allpats(:,pat),segMWCTLP(:,pat),'type','spearman');
    fit2 = polyfit(segWork_allpats(:,pat),segMWCTLP(:,pat),1);
    text(1.5,-3.2,['r^2 = ',num2str(r2.^2,'%.2f')])
    text(1.5,-4,['y = ',num2str(fit2(1),'%.2f'),'x + ',num2str(fit2(2),'%.2f')])
    axis([-0.5 4 -5 4],'square')
    ylabel({'Static Laplace';'Stress-Strain Area via CT'})
    xlabel('Stress-Strain Area via PSM')
    
    subplot(1,3,3); hold all
    plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
    plot(segWork_allpats(:,pat),segMWCTLPTV(:,pat),'b.');
    [r3,p] = corr(segWork_allpats(:,pat),segMWCTLPTV(:,pat),'type','spearman');
    fit3 = polyfit(segWork_allpats(:,pat),segMWCTLPTV(:,pat),1);
    text(1.5,-3.2,['r^2 = ',num2str(r3.^2,'%.2f')])
    text(1.5,-4,['y = ',num2str(fit3(1),'%.2f'),'x + ',num2str(fit3(2),'%.2f')])
    axis([-0.5 4 -5 4],'square')
    ylabel({'Time-Varying Laplace'; 'Stress-Strain Area via CT'})
    xlabel('Stress-Strain Area via PSM')
    sgtitle(['LV Regional MW Agreement for Patient ',num2str(pat)])
    set(findall(gcf,'-property','FontSize'),'FontSize',22)
    set(findall(gcf,'-property','LineWidth'),'LineWidth',4)
    set(findall(gcf,'-property','MarkerSize'),'MarkerSize',38)

    r_allpats(pat,:) = [r1,r2,r3];
    m_allpats(pat,:) = [fit1(1),fit2(1),fit3(1)];
end

r_stats = [median(r_allpats);prctile(r_allpats,25);prctile(r_allpats,75)];
m_stats = [median(m_allpats);prctile(m_allpats,25);prctile(m_allpats,75)];

%%
patcolors = turbo(8);
patcolors2 = [patcolors(2:8,:);patcolors(1,:)];
figure;
for pat = 1:8
    subplot(2,4,pat); hold all; 
    plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
    plot(segWork_allpats(:,pat),segMWCTLPTV(:,pat),'.','Color',patcolors2(pat,:));
    [r3,p] = corr(segWork_allpats(:,pat),segMWCTLPTV(:,pat),'type','spearman');
    fit3 = polyfit(segWork_allpats(:,pat),segMWCTLPTV(:,pat),1);
    text(0.8,-3.2,['r = ',num2str(r3,'%.2f')])
    text(0.8,-4,['y = ',num2str(fit3(1),'%.2f'),'x + ',num2str(fit3(2),'%.2f')])
    axis([-0.5 4 -5 4],'square')
    ylabel({'Time-Varying Laplace'; 'Stress-Strain Area via CT'})
    xlabel('Stress-Strain Area via PSM')
    sgtitle('LV Regional MW Agreement')
    title(['Patient ',num2str(pat)])
    set(findall(gcf,'-property','FontSize'),'FontSize',22)
    set(findall(gcf,'-property','LineWidth'),'LineWidth',4)
    set(findall(gcf,'-property','MarkerSize'),'MarkerSize',38)
    set(findall(gcf,'-property','Box'),'Box','on')
end

%% Segment break down
figure; 
for region = 1:17
    subplot(4,5,region); hold all;
    plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
    plot(segWork_allpats,segMWCTLPTV,'.','Color',[0.7 0.7 0.7]);
    plot(segWork_allpats(region,:),segMWCTLPTV(region,:),'r.');
    ylabel({'Time-Varying Laplace'; 'Stress-Strain Area via CT'})
    xlabel('Stress-Strain Area via PSM')
    title(['Segment ',num2str(region)])
    axis([-0.5 4 -5 4])
end
set(findall(gcf,'-property','FontSize'),'FontSize',18)
set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',25)

%% Maybe test and see if you plot all the segmental ellipsoid fits in the same plot, do they all go in the right direction?
for pat = 5
    time = 1;
    foldpath = ['/Users/amandacraine/Documents/ContijochLab/BiVFunctionModels/BiV1-8/BiV',num2str(pat)];
    modellist = dir([foldpath,'/*.obj']);
    model = readObj([foldpath,'/',modellist(time).name]); %End-diastole
    framepts = model.v; %individuals vertex pts
    triangulation = model.f.v; %faces
    thickness_data = readmatrix([datapath, 'BiV',num2str(pat),'/WallThicknessWholeLVTimeVar.txt']);

    elemslist = dir([foldpath,'/*.txt']);
    elems = readtable([elemslist(1).folder,'/',elemslist(time).name]);
    ElemList = table2array(elems(:,1));

    Elem2Tri = zeros(size(triangulation));
    for i = 1:size(triangulation,1)
        for j = 1:size(triangulation,2)
            Elem2Tri(i,j) = ElemList(triangulation(i,j));
        end
    end

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

    figure; %set(gcf,'Position',[100 100 1500 1500]) %hold all
    for region = 1:17
        x = framepts(triangulation(TriSegElems == region,:),3);
        y = framepts(triangulation(TriSegElems == region,:),2);
        z = -framepts(triangulation(TriSegElems == region,:),1);

        [center, radii(:,region), evecs, v, chi2] = ellipsoid_fit([x,y,z], '0xy');


        regionLoc = [4,1;
            1,2;
            1,6;
            4,7;
            7,6;
            7,2;
            4,2;
            2,3;
            2,5;
            4,6;
            6,5;
            6,3;
            4,3;
            3,4;
            4,5;
            5,4;
            4,4];

        subplot(7,7,(regionLoc(region,2)-1)*7+regionLoc(region,1));

        plot3( x, y, z, '.r' );
        %hold on;

        %draw fit
        mind = min( [ x y z ] );
        maxd = max( [ x y z ] );
        nsteps = 50;
        step = ( maxd - mind ) / nsteps;
        [ x, y, z ] = meshgrid( linspace( mind(1) - step(1), maxd(1) + step(1), nsteps ), linspace( mind(2) - step(2), maxd(2) + step(2), nsteps ), linspace( mind(3) - step(3), maxd(3) + step(3), nsteps ) );

        Ellipsoid = v(1) *x.*x +   v(2) * y.*y + v(3) * z.*z + ...
            2*v(4) *x.*y + 2*v(5)*x.*z + 2*v(6) * y.*z + ...
            2*v(7) *x    + 2*v(8)*y    + 2*v(9) * z;
        p = patch(isosurface( x, y, z, Ellipsoid, -v(10)));
        hold off;
        set( p, 'FaceColor', 'g', 'EdgeColor', 'none' );
        view( -70, 40 );
        axis vis3d equal;
        camlight;
        lighting phong;
        hold on;
    % 
    % radius = (radii(1,region).^2).*(radii(2,region).^2).*(radii(3,region).^2).*(x.^2./radii(1,region).^4 + y.^2./radii(2,region).^4 + z.^2./radii(3,region).^4).^2;
    % r_princ(region,pat) = mean(radius);
    end

    r_eff = radii(1,:).*radii(2,:)./(radii(1,:) + radii(2,:));

    
    %r_eff = radii(1,:).*radii(2,:).*radii(3,:)./(radii(2,:).*radii(3,:) + radii(1,:).*radii(3,:) + radii(1,:).*radii(2,:));
    %LPtest = repmat(LVP_cath{pat},17,1).*r_eff'./2./thickness_data(:,1);
    MWLPtest(:,pat) = calculateSegmentalWorkEstimates(RS_CT{pat},LVP_cath{pat},r_princ(:,pat),thickness_data(:,1),elementIndices,ElemList,triangulation,17,t{pat});
end
% load all_seg_work_all_pats.mat
% load segMWCT_effrad_allpats.mat


figure;
subplot(1,2,1); set(gca,'colororder',turbo(8)); hold all
plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
plot(segWork_allpats,segMWCTLP,'.');
[r,p] = corr(reshape(segWork_allpats,numel(segWork_allpats),1),reshape(segMWCTLP,numel(segMWCTLP),1),'type','spearman');
fitvals = fit(reshape(segWork_allpats,numel(segWork_allpats),1),reshape(segMWCTLP,numel(segMWCTLP),1),'poly1');
text(1.5,-3.2,['r = ',num2str(r,'%.2f')])
text(1.5,-4,['y = ',num2str(fitvals.p1,'%.2f'),'x + ',num2str(fitvals.p1,'%.2f')])
%axis([-0.5 4 -5 4],'square')
ylabel({'Static Laplace';'Stress-Strain Area via CT'})
xlabel('Stress-Strain Area via PSM')

subplot(1,2,2); set(gca,'colororder',turbo(8)); hold all
plot([-5 4],[-5 4],'--','Color',[0.9 0.9 0.9])
plot(segWork_allpats,MWLPtest,'.');
[r,p] = corr(reshape(segWork_allpats,numel(segWork_allpats),1),reshape(MWLPtest,numel(MWLPtest),1),'type','spearman');
fitvals = fit(reshape(segWork_allpats,numel(segWork_allpats),1),reshape(MWLPtest,numel(MWLPtest),1),'poly1');
text(1.5,-3.2,['r = ',num2str(r,'%.2f')])
text(1.5,-4,['y = ',num2str(fitvals.p1,'%.2f'),'x + ',num2str(fitvals.p2,'%.2f')])
%plot(fitvals,reshape(segWork_allpats,numel(segWork_allpats),1),reshape(MWLPtest,numel(MWLPtest),1))
%axis([-0.5 4 -5 4],'square')
ylabel({'Static Laplace';'Stress-Strain Area via CT'})
xlabel('Stress-Strain Area via PSM')




