function [NegFraction] = calculateNegWorkFraction(volume,work)
%given local volume and myocardial work information, we find the locations
%where work < 0. Then the negfraction considers the fraction of the LV that
%is negative given the volume of the negative locations.

% lvVolume			= detj{pat}(gauss{1});
% lvWork				= work{pat}(gauss{1});
NegLocation		= find(work < 0.0);
NegVolume			= sum(volume(NegLocation));
NegFraction	= NegVolume/sum(volume);



end