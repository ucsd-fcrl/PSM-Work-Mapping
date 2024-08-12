clc;
%cd('/Volumes/MainShare/projects/PSM/PSM/')
numPatients = 8;
%patOrder = [6 3 5 1 8 2 4 7];
patOrder = 1:8;

lvNegFraction = cell(1,numsol);
rvNegFraction = cell(1,numsol);
sepNegFraction = cell(1,numsol);
lvfwNegFraction = cell(1,numsol);
lvFibNegFraction = cell(1,numsol);
fibNegFraction = cell(1,numsol);
intersectFraction = cell(1,numsol);
negFraction = cell(1,numsol);
totalWork = cell(1,numsol);
medianWork = cell(1,numsol);
medianNegWork = cell(1,numsol);

stresspfElem = cell(numsol,1);
stresspsElem = cell(numsol,1);
stresspcElem = cell(numsol,1);
stretchpfElem = cell(numsol,1);
stretchpsElem = cell(numsol,1);
stretchpcElem = cell(numsol,1);
workElem = cell(numsol,1);

for i=1:numsol
	medianWork{i} = median(workVol{i});
	medianNegWork{i} = median(workNegVol{i});
	totalVolume = sum(detj{i});
	scarVolume{i} = sum(detj{i}.*infarctGauss{i})/(sum(detj{i}(gauss{1})));

	lv90Work{i} = prctile(workVol{i}(gauss{1}),90);
	totalNegWork{i} = sum(workNeg{i}.*detj{i});
	totalNegWorkDensity{i} = totalNegWork{i}/totalVolume;
	totalWork{i} = sum(workVol{i});
	totalWorkDensity{i} = totalWork{i}/totalVolume;
	
	negFraction{i}		= length(find(work{i} < 0.0))/length(work{i});
	lvNegFraction2{i}	= length(find(work{i}(gauss{1}) < 0.0))/length(gauss{1});
	sepNegFraction2{i}	= length(find(work{i}(gauss{3}) < 0.0))/length(gauss{3});
	rvNegFraction2{i}	= length(find(work{i}(gauss{2}) < 0.0))/length(gauss{2});
	lvfwNegFraction2{i}	= length(find(work{i}(gauss{4}) < 0.0))/length(gauss{4});

	myo75Work			= prctile(work{i}, 75);
	high75Location		= find(work{i} > myo75Work);
	high75Volume		= sum(detj{i}(high75Location));
	high75Fraction{i}	= high75Volume/totalVolume;
	
	
% 	negLocation			= find(work{i} < 0.0);
% 	negVolume		= sum(detj(negLocation));
		
	lvVolume			= detj{i}(gauss{1});
	lvWork				= work{i}(gauss{1});
	lv75Work			= prctile(lvWork,75);
	lvNegLocation		= find(lvWork < 0.0);
	lv75Location		= find(lvWork > lv75Work);
	lvNegVolume			= sum(lvVolume(lvNegLocation));
	lv75Volume			= sum(lvVolume(lv75Location));
	lvNegFraction{i}	= lvNegVolume/sum(lvVolume);
	lv75Fraction{i}		= lv75Volume/sum(lvVolume);
	
	lvNegLocation		= find(lvWork < 0.0);
	lvNegVolume			= sum(lvVolume(lvNegLocation));
	lvNegFraction{i}	= lvNegVolume/sum(lvVolume);
	
	rvVolume			= detj{i}(gauss{2});
	rvWork				= work{i}(gauss{2});
	rvNegLocation		= find(rvWork < 0.0);
	rvNegVolume			= sum(rvVolume(rvNegLocation));
	rvNegFraction{i}	= rvNegVolume/sum(rvVolume);
	
	sepVolume			= detj{i}(gauss{3});
	sepWork				= work{i}(gauss{3});
	sepNegLocation		= find(sepWork < 0.0);
	sepNegVolume		= sum(sepVolume(sepNegLocation));
	sepNegFraction{i}	= sepNegVolume/sum(sepVolume);
	
	lvfwVolume			= detj{i}(gauss{4});
	lvfwWork			= work{i}(gauss{4});
	lvfwNegLocation		= find(lvfwWork < 0.0);
	lvfwNegVolume		= sum(lvfwVolume(lvfwNegLocation));
	lvfwNegFraction{i}	= lvfwNegVolume/sum(lvfwVolume);
		
	lvFibNegFraction{i} = length(find(workF{i}(gauss{1}) < 0.0))/length(gauss{3});
 	fibNegFraction{i}	= length(find(workF{i} < 0.0))/length(work{i});
end
negFraction
lvNegFraction
rvNegFraction
sepNegFraction

meanWorkt = meanWork';
stdWorkt = stdWork';
covWork = (stdWork./meanWork)';
meanNegWorkt = meanNegWork';
stdNegWorkt = stdNegWork';
covNegWork = (stdNegWork./meanNegWork)';

meanWorkt(1,:)
stdWorkt(1,:)
covWork(1,:)

CRTModels = true; %true
if (CRTModels)
	p = numPatients;
	lvNegFractionValues = cell2mat(lvNegFraction);
	lvNegFractionImpFrac = (lvNegFractionValues(1:p) - lvNegFractionValues(p+1:2*p))./lvNegFractionValues(1:p);
	lvNegFractionImp = (lvNegFractionValues(1:p) - lvNegFractionValues(p+1:2*p));

	rvNegFractionValues = cell2mat(rvNegFraction);
	rvNegFractionImpFrac = (rvNegFractionValues(1:p) - rvNegFractionValues(p+1:2*p))./rvNegFractionValues(1:p);
	rvNegFractionImp = (rvNegFractionValues(1:p) - rvNegFractionValues(p+1:2*p));

	sepNegFractionValues = cell2mat(sepNegFraction);
	sepNegFractionImpFrac = (sepNegFractionValues(1:p) - sepNegFractionValues(p+1:2*p))./sepNegFractionValues(1:p);
	sepNegFractionImp = (sepNegFractionValues(1:p) - sepNegFractionValues(p+1:2*p));

	lvfwNegFractionValues = cell2mat(lvfwNegFraction);
	lvfwNegFractionImpFrac = (lvfwNegFractionValues(1:p) - lvfwNegFractionValues(p+1:2*p))./lvfwNegFractionValues(1:p);
	lvfwNegFractionImp = (lvfwNegFractionValues(1:p) - lvfwNegFractionValues(p+1:2*p));

	negFractionValues = cell2mat(negFraction);
	negFractionImpFrac = (negFractionValues(1:p) - negFractionValues(p+1:2*p))./negFractionValues(1:p);
	negFractionImp = (negFractionValues(1:p) - negFractionValues(p+1:2*p));

	covNegWorkImp = covNegWork(:,1:p) - covNegWork(:,p+1:2*p);
	covNegWorkImpFrac = covNegWorkImp./covNegWork(:,1:p);

	covWorkImp = covWork(:,1:p) - covWork(:,p+1:2*p);
	covWorkImpFrac = covWorkImp./covWork(:,1:p);

	totalNegWorkValues = cell2mat(totalNegWork);
	totalNegWorkImp = totalNegWorkValues(1:p) - totalNegWorkValues(p+1:2*p);
	totalNegWorkImpFrac = totalNegWorkImp./totalNegWorkValues(1:p);
	
	if (compISF)
		ISFOuterImp = ISFOuter(1:p) - ISFOuter(p+1:2*p);
		ISFOuterImpFrac = ISFOuterImp./ISFOuter(1:p);
		ISFInnerImp = ISFInner(1:p) - ISFInner(p+1:2*p);
		ISFInnerImpFrac = ISFInnerImp./ISFInner(1:p);
		ISFImp = ISF(1:p) - ISF(p+1:2*p);
		ISFImpFrac = ISFImp./ISF(1:p);
	end	
end

findSelfIntersection = false;
if (findSelfIntersection)
	lvrvsep = 1; % 1 = lv, 2 = rv ,3 = sepv
	totalGaussPoints = numGaussPoints*numElements;
	for i=1:numsol
		numIntersectingGaussPoints = 0;
		for gaussPtNum = 1:length(gauss{lvrvsep})
			currentGaussPoint = gauss{lvrvsep}(gaussPtNum);
			[xInt,yInt,segments]= selfintersect(log(stretchpf{i}(currentGaussPoint,:)'),stresspf{i}(currentGaussPoint,:));
			numLoops = 0;
			for seg = 1:size(segments,1)
				currentSeg = segments(seg,:);
				if (currentSeg(2) < (length(t{i}) - 20))
					numLoops = numLoops+1;
				end
			end
			if (numLoops > 0)
				numIntersectingGaussPoints = numIntersectingGaussPoints+1;
			end
		end
		intersectFraction{i} = numIntersectingGaussPoints/length(gauss{lvrvsep});
	end
end
