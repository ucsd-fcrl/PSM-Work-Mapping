function [RS_CT,patch_areas] = calculate_RSCT(foldpath,modellist,triangulation,areas)
% Calculate FAC for triangulations aka SQUEEZ
for i = 1:length(modellist)
    model = readObj([foldpath,'/',modellist(i).name]);
    framepts = model.v; %individuals vertex pts

    A = framepts(triangulation(:,1),:);
    B = framepts(triangulation(:,2),:);
    C = framepts(triangulation(:,3),:);
    AB = B - A;
    lengthAB = sqrt(sum(AB.^2,2));
    AC = C - A;
    lengthAC = sqrt(sum(AC.^2,2));
    BC = B - C;
    lengthBC = sqrt(sum(BC.^2,2));
    s = (lengthAB + lengthBC + lengthAC)/2;
    patch_areas(:,i) = sqrt(s.*(s-lengthAB).*(s-lengthAC).*(s-lengthBC));
end

normalized_areas = sqrt(patch_areas ./ patch_areas(:,1));
RS_CT = normalized_areas - 1;

end