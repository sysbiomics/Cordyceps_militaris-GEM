function [growthProfile] = testPOPCORN(model,biomassRxn,targetRxn,carbonSource,nitrogenSource,nCarbon,nNitrogen)
% POPCORN was developed for determining the optimal PrOPortion of CarbOn and 
% nitRogeN for growth and target production 
% 
% USAGE:
%    [growthProfile] = POPCORN(model,,biomassRxn,targetRxn,carbonSource,nitrogenSource,nCarbon,nNitrogen)
%
% INPUTS:
%    model:            RAVEN model structure
%    biomassRxn:       reaction ID of the biomass production or growth reaction
%    targetRxn:        reaction ID of target production reaction
%    carbonSource:     reaction ID of uptake reaction of carbon source
%    nitrogenSource:   reaction ID of uptake reaction of nitrogen source
%    nCarbon:          number of carbon atom(s) of carbon source
%    nNitrogen:        number of nitrogen atom(s) of nitrogen source
%
% OUTPUTS:
%    growthProfile:    a structure of: 1) cell of C:N ratio series, named CNratios
%                                      2) cell of the maximum growth rate (/h), named umax
%                                      3) cell of the maximum target production (mmol/gDW/h), named targetP
%                                      4) cell of the yield of target on
%                                      biomass (mmol/gDW), named targetY
%                                      5) cell of the upper bound used for
%                                      the uptake reaction of carbon
%                                      source, named carbonLB
%                                      6) cell of the upper bound used for
%                                      the uptake reaction of nitrogen
%                                      source, named nitrogenLB
%
%
% Written by Nachon Raethong, 02-AUG-2019
% Updated by Nachon Raethong, 08-AUG-2019
%
%

% generate a series of upper bounds used for the uptake reaction of nitrogen source
variedN = [0.01:0.2:19.99]; % 100 values

% create a structure and cells for the results
growthProfile = struct();
growthProfile.CNratios = cell(numel(variedN),1);
growthProfile.umax = cell(numel(variedN),1);
growthProfile.targetP = cell(numel(variedN),1);
growthProfile.targetY = cell(numel(variedN),1);
growthProfile.carbonUP = cell(numel(variedN),1);
growthProfile.nitrogenUP = cell(numel(variedN),1);

%Find out the umax, target production and yield in different C:N ratios iteratively
model = setParam(model,'obj',biomassRxn,1);
for i = 1:numel(variedN)
    % generate an upper bound used for the uptake reaction of carbon source
    variedC = 20 - variedN(i);
    % define the upper bounds for the uptake reactions of carbon and
    % nitrogen sources
    model = setParam(model,'ub',carbonSource,variedC);
    model = setParam(model,'ub',nitrogenSource,variedN(i));
    %Find out the umax 
    sol = solveLP(model,1); 
        
    %Find out the target production with 90% of umax
    model2 = setParam(model,'lb',biomassRxn,(sol.f*-0.9));
    model2 = setParam(model2,'obj',targetRxn,1);
    sol2 = solveLP(model2,1);
    
    %Generating output
    growthProfile.CNratios{i} = ((variedC*nCarbon)/(variedN(i)*nNitrogen));
    growthProfile.umax{i} = (sol.f*-1);
    growthProfile.targetP{i} = (sol2.f*-1);
    growthProfile.targetY{i} = ((sol2.f*-1)/(sol.f*-0.9));
    growthProfile.carbonUP{i} = (variedC);
    growthProfile.nitrogenUP{i} = (variedN(i));
    
end
%Show the output on screen
struct2table(growthProfile)

%Draw the resulting scatter plots of umax and targetP againt CNratio series
x = cell2mat(growthProfile.CNratios);
y1 = cell2mat(growthProfile.umax);
y2 = cell2mat(growthProfile.targetP);
createfigure(x,y1,15,[1 0 1],y2,[0 0 0]);

end
