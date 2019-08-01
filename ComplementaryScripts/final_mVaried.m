%% Growth capability of C. militaris underlying nutritional deprivation
% Purpose: to determine growth capability of C. militaris underlying nutritional deprivation using iNR1320
% 
% Written by Nachon Raethong, 31-JAN-2019
% Updated by Nachon Raethong, 24-JULY-2019

%% WORKSPACE
cd 'C:\Users\Dell\Documents\GitHub\Cordyceps_militaris-GEM';


load('C:\Users\Dell\Documents\GitHub\Cordyceps_militaris-GEM\ModelFiles\mat\model.mat');
modelBigger = model;
sol = solveLP(modelBigger,1);
printFluxes(modelBigger,sol.x);

modelRaven = modelBigger;
uptakeRxns = modelRaven.rxns(strncmp(modelRaven.rxnNames,'uptake ',7));
uptakeRxnNames = modelRaven.rxnNames(strncmp(modelRaven.rxnNames,'uptake ',7));
requiredRxns = {'o2IN' 'piIN' 'slfIN'};

modelRaven = setParam(modelRaven,'ub',uptakeRxns,0); % block all uptake
modelRaven = setParam(modelRaven,'ub',requiredRxns,1000); % allow only required metabolites
modelRaven = setParam(modelRaven,'lb',{'bmOUT' 'cordycepinOUT'},0);
modelRaven = setParam(modelRaven,'ub',{'bmOUT' 'cordycepinOUT'},1000);
modelRaven = setParam(modelRaven,'obj',{'bmOUT'},1);

mVersatile_I = cell(numel(uptakeRxns),1);
solVersatile_I= cell(numel(uptakeRxns),1);
model_mVersatile_I = modelRaven; % allow only required metabolites 
for i = 1:numel(uptakeRxns)
    model = model_mVersatile_I;
    model_i = setParam(model,'ub',uptakeRxns{i},2);
    sol = solveLP(model_i,1);
    mVersatile_I{i} = (sol.f*-1);
    solVersatile_I{i} = sol;
end

solVersatile_II= cell(numel(uptakeRxns),1);
mVersatile_II = cell(numel(uptakeRxns),1);
model_mVersatile_II = setParam(modelRaven,'ub','glcIN',1); % allow only required metabolites and glucose
for i = 1:numel(uptakeRxns)
    model = model_mVersatile_II;
    model_i = setParam(model,'ub',uptakeRxns{i},1);
    sol = solveLP(model_i,1);
    mVersatile_II{i} = (sol.f*-1);
    solVersatile_II{i} = sol;
end


solVersatile_III = cell(numel(uptakeRxns),1);
mVersatile_III = cell(numel(uptakeRxns),1);
model_mVersatile_III = setParam(modelRaven,'ub','nh3IN',1); % allow only required metabolites and ammonia
for i = 1:numel(uptakeRxns)
    model = model_mVersatile_III;
    model_i = setParam(model,'ub',uptakeRxns{i},1);
    sol = solveLP(model_i,1);
    mVersatile_III{i} = (sol.f*-1);
    solVersatile_III{i} = sol;
end

variedResult = table(uptakeRxns,uptakeRxnNames,mVersatile_I,mVersatile_II,mVersatile_III);

% The result was used to generate a heatmap representing the prediction growth rates for 46 individual 
% nutritional substrates governed under nutrient deprivation, nitrogen source deprivation and carbon source 
% deprivation (Fig 3).


