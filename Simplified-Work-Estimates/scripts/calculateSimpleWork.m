%this script calculates the four simplified myocaridal work estmiates for
%all patches, not just the segments. 
% segmental measurements like radius and wall thickness are applied to all
% points that fall into one segment

%we're calculating the laplace work measurements and work with generic
%pressure waveforms.
clear; clc

% NOTE: make sure your home path is the Simplified Work Estimates directory
homepath = '/Users/amandacraine/Documents/ContijochLab/repos/PSM-Work-Mapping/Simplified-Work-Estimates';
addpath([homepath,'/data/'])
addpath([homepath,'/scripts'])
addpath([homepath,'/LV Geometric Models'])
cd(homepath)
savepath = [homepath '/data-collected/'];
load RS_CT.mat
load dataPSMLBBBCRT.mat

%% Laplace wall stress-strain area
load segMWCT_effrad_allpats.mat
load Laplace_measurements.mat
warning off
%get segmental lapalce calculations
for pat = 1:8
    disp(['Analyzing patient ',num2str(pat)])
    foldpath = [homepath '/LV Geometric Models/BiV1-8/BiV',num2str(pat)];
    for time = 1:length(t{pat})
        %load in thickness data
        modellist = dir([foldpath,'/*.obj']);
        model = readObj([foldpath,'/',modellist(time).name]); %End-diastole
        framepts = model.v; %individuals vertex pts
        triangulation = model.f.v; %faces
        %thickness_data = readmatrix([datapath, 'BiV',num2str(pat),'/WallThicknessWholeLVTimeVar.txt']);
    
        elemslist = dir([foldpath,'/*.txt']);
        elems = readtable([elemslist(1).folder,'/',elemslist(time).name]);
        ElemList = table2array(elems(:,1));
    end

    MWCTLP(:,pat) = calculatePatchWorkEstimates(RS_CT{pat},LVP_cath{pat},r_eff_allpats{pat}(:,1),thickness_data{pat}(:,1),elementIndices,ElemList,triangulation,17,t{pat});
    MWCTLPTV(:,pat) = calculatePatchWorkEstimates(RS_CT{pat},LVP_cath{pat},r_eff_allpats{pat},thickness_data{pat},elementIndices,ElemList,triangulation,17,t{pat});
    

end
save([savepath 'Laplace_work_all_patches.mat'],'MWCTLP','MWCTLPTV')
%apply to each patch in their respective segment

%% generic pressure waveforms

load genericLVPs.mat

for pat = 1:8
    disp(['Analyzing patient ',num2str(pat)])
    foldpath = [homepath '/LV Geometric Models/BiV1-8/BiV',num2str(pat)];
    for time = 1:length(t{pat})
        %load in thickness data
        modellist = dir([foldpath,'/*.obj']);
        model = readObj([foldpath,'/',modellist(time).name]); %End-diastole
        framepts = model.v; %individuals vertex pts
        triangulation = model.f.v; %faces
        %thickness_data = readmatrix([datapath, 'BiV',num2str(pat),'/WallThicknessWholeLVTimeVar.txt']);
    
        elemslist = dir([foldpath,'/*.txt']);
        elems = readtable([elemslist(1).folder,'/',elemslist(time).name]);
        ElemList = table2array(elems(:,1));
    end

    MWCTgenericLVP(:,pat) = calculatePatchWorkEstimates(RS_CT{pat},genericLVP{pat},[],[],elementIndices,ElemList,triangulation,17,t{pat});
    MWCTgenericLVP_patES(:,pat) = calculatePatchWorkEstimates(RS_CT{pat},genericLVP_ESP{pat},[],[],elementIndices,ElemList,triangulation,17,t{pat});
end
save([savepath 'GenericLVP_work_all_patches.mat'],'MWCTgenericLVP','MWCTgenericLVP_patES')




