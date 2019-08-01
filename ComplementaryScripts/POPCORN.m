function growthProfile = POPCORN(model)
% This is the constain-based modeling function, POPCORN used to determine the optimal PrOPortion of CarbOn and 
% nitRogeN for growth and cordycepin production in C. militaris 
%
% INPUT:
%		model = model (in this case, it must be iNR1320)
% OUTPUT:
%		growthProfile = the table of CNratios and their umax, cordycepin production, glucose uptake rate and 
% ammonia uptake rate
% 
% Usage:
% 			growthProfile = POPCORN(model)
%
% Written by Nachon Raethong, 02-AUG-2019
%
%
sugars = {'glcIN' 'fruIN' 'xylIN' 'arabIN' 'sucIN'};
model = setParam(model,'ub',sugars,0);
model = setParam(model,'ub','cordycepinOUT',1000);
model = setParam(model,'lb','cordycepinOUT',0);
model = setParam(model,'obj','bmOUT',1);

variedN = [0.05:0.04:1.05];
CNratios = cell(numel(variedN),1);
umax = cell(numel(variedN),1);
cordycepin = cell(numel(variedN),1);
glucose = cell(numel(variedN),1);
ammonia = cell(numel(variedN),1);

for i = 1:numel(variedN)
    variedC = 1.1 - variedN(i);
    model = setParam(model,'ub','glcIN',variedC);
    model = setParam(model,'ub','nh3IN',variedN(i));
    sol = solveLP(model,1);
    CNratios{i} = num2str((variedC*6)/(variedN(i)*1));
    umax{i} = num2str((sol.f*-1)/1.1);
    model2 = setParam(model,'lb','bmOUT',((sol.f*-1)/1.125));
    model2 = setParam(model2,'obj','cordycepinOUT',1);
    sol2 = solveLP(model2,1);
    cordycepin{i} = num2str(sol2.f*-1);
    glucose{i} = variedC;
    ammonia{i} = variedN(i);

end

growthProfile = table(CNratios,umax,cordycepin,glucose,ammonia);
end