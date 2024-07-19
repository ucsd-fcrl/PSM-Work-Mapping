function [NegWorkFraction] = getNegWorkFractions(work,areas)

    workNegLocation   = find(work < 0);
    workNegVolume     = sum(areas(workNegLocation));
    NegWorkFraction = workNegVolume./sum(areas);

end