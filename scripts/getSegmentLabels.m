function [pointsSeg,TriSegElems] = getSegmentLabels(elementIndices,pat,endo_flag)
%endo_flag is to indicate if we're analyzing only the LV endocardium (1) or the
%whole LV (0)

if endo_flag == 1
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
else
    segmentElems{2} = [1,2,3,4];
    segmentElems{1} = [5,7,69,71];
    segmentElems{6} = [6,8,70,72];
    segmentElems{5} = [9,11,73,75];
    segmentElems{4} = [10,12,74,76];
    segmentElems{3} = [13,14,15,16];

    segmentElems{8} = [18,20,22,24];
    segmentElems{7} = [26,30,90,94];
    segmentElems{12} = [28,32,92,96];
    segmentElems{11} = [34,38,98,102];
    segmentElems{10} = [36,40,100,104];
    segmentElems{9} = [42,44,46,48];

    segmentElems{14} = [17,19,41,43,21,23,45,47];
    segmentElems{15} = [35,39,99,103];
    segmentElems{16} = [27,33,31,37,91,97,95,101];
    segmentElems{13} = [25,29,89,93];
    segmentElems{17} = [49,50,53,54,57,58,61,62,51,52,55,56,59,60,63,64,...
        113,114,117,118,121,122,125,126,115,116,119,120,123,124,127,128];
end


pointsSeg = zeros(length(elementIndices),1);
for i = 1:length(elementIndices)
    for k = 1:length(segmentElems)
        if any((segmentElems{k} == elementIndices(i)))
            pointsSeg(i) = k;
        end
    end
end


%apply mean segmental volume to all points
foldpath = ['/Users/amandacraine/Documents/ContijochLab/BiVFunctionModels/BiV1-8/BiV',num2str(pat)];
modellist = dir([foldpath,'/*.obj']);
model = readObj([foldpath,'/',modellist(1).name]);
triangulation = model.f.v; %faces

elemslist = dir([foldpath,'/*.txt']);
elems = readtable([elemslist(1).folder,'/',elemslist(1).name]);
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


end