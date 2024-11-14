%this script computes the area of each patch
clear; clc
% NOTE: make sure your home path is the Simplified Work Estimates directory
homepath = '/Users/amandacraine/Documents/ContijochLab/repos/PSM-Work-Mapping/Simplified-Work-Estimates';
addpath([homepath,'/data/'])
addpath([homepath,'/scripts'])
addpath([homepath,'/LV Geometric Models'])
cd(homepath)
load dataPSMLBBBCRT.mat
savepath = [homepath '/data-collected/'];
PatchAreas = cell(8,1);
for pat = 1:8
    foldpath = [homepath '/LV Geometric Models/BiV1-8/BiV',num2str(pat)];
    modellist = dir([foldpath,'/*.obj']);
    model = readObj([foldpath,'/',modellist(1).name]);
    triangulation = model.f.v; %faces
    areas = zeros(size(triangulation,1),length(modellist));
    [~,patch_areas] = calculate_RSCT(foldpath,modellist,triangulation,areas);
    PatchAreas{pat} = patch_areas;

    [pointsSeg_endo,~] = getSegmentLabels(elementIndices,pat,1);
    endo_idx = find(pointsSeg_endo > 0);
    endo_vols = detj{pat}(endo_idx);
    [sum(detj{pat}(gauss{1})) sum(patch_areas(:,1))]
end



save([savepath 'PatchAreas.mat'],'PatchAreas')