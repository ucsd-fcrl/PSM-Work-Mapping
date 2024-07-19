function [work_estimate,seg_strain,seg_stress,meanwork] = calculateSegmentalWorkEstimates(strain_estimate,stress_estimate,radius,thickness,elementIndices,ElemList,triangulation,numsegs,time)
%%Purpose: calculate work estimates for each segment given a stress and
%%strain input
    
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

    %%Compare MW
    
    pointsSeg = zeros(length(elementIndices),1);
    for i = 1:length(elementIndices)
        for k = 1:length(segmentElems)
            if any((segmentElems{k} == elementIndices(i)))
                pointsSeg(i) = k;
            end
        end
    end




seg_radius = zeros(numsegs,length(time));
seg_strain = zeros(numsegs,length(time));
seg_stress = zeros(numsegs,length(time));
work_estimate = zeros(numsegs,1);

for region = 1:numsegs
    %to get the segmental PSM work only
    if isempty(strain_estimate)
        work_estimate(region,:) = mean(stress_estimate(pointsSeg == region,:));

        % for my MW estimates:
        % 1)  for pressure-strain estimate, does not include shape info
    elseif isempty(radius) && isempty(thickness)
        seg_strain(region,:) = mean(strain_estimate(TriSegElems == region,:));
        %seg_stress(region,:) = repmat(stress_estimate,17,1);
        work_estimate(region,:) = PolyAreaSigned(seg_strain(region,:),stress_estimate);

        for l = 1:1024
            work_by_patch(l) = PolyAreaSigned(strain_estimate(l,:),stress_estimate);
        end
        meanwork = mean(work_by_patch);
        % 2)  for shape-informed wall stress estimates
    elseif  size(radius,1) > numsegs
        %for the principal radius measurements, need to find average
        %principal radius for each segment
        seg_strain(region,:) = mean(strain_estimate(TriSegElems == region,:));
        seg_radius(region,:) = mean(radius{pat}(TriSegElems == region,:));
        seg_stress(region,:) = stress_estimate.*seg_radius(region,:)./(2.*thickness(region,:));
        % seg_stress(region,:) = stress_estimate.*seg_radius(region,:)./...
        %     (thickness(region,:).*(2 + (thickness(region,:)./seg_radius(region,:))));
        work_estimate(region,:) = PolyAreaSigned(seg_strain(region,:),seg_stress(region,:));
        meanwork = [];%mean(work_by_patch);
        %3)  for the radius measurements acquired by the model
    elseif size(radius,1) <= numsegs
        seg_strain(region,:) = mean(strain_estimate(TriSegElems == region,:));
        seg_stress(region,:) = stress_estimate.*radius(region,:)./(2.*thickness(region,:));
        % seg_stress(region,:) = stress_estimate.*radius(region,:)./...
        %     (thickness(region,:).*(2 + (thickness(region,:)./radius(region,:))));
        work_estimate(region,:) = PolyAreaSigned(seg_strain(region,:),seg_stress(region,:));
        meanwork = [];%mean(work_by_patch);

    end

    %work_estimate(region,:) = PolyAreaSigned(seg_strain(region,:),seg_stress(region,:));
    
end

end