cd 'C:\Users\Dell\Documents\GitHub\cordyceps_model';
load('C:\Users\Dell\Documents\GitHub\cordyceps_model\ModelFiles\mat\model.mat');
modelBigger = model;
sol = solveLP(modelBigger,1);
printFluxes(modelBigger,sol.x);


modelRaven = modelBigger;
%uptakeRxns = {'glcIN' 'nh3IN' 'sucIN' 'glnIN' 'chitIN'};
uptakeRxns = modelRaven.rxns(strncmp(modelRaven.rxnNames,'uptake ',7));
uptakeRxnNames = modelRaven.rxnNames(strncmp(modelRaven.rxnNames,'uptake ',7));
requiredRxns = {'o2IN' 'piIN' 'slfIN'};

modelRaven = setParam(modelRaven,'ub',uptakeRxns,0); % block all uptake
%modelRaven = setParam(modelRaven,'lb',uptakeRxns,-1000); % block all uptake
modelRaven = setParam(modelRaven,'ub',requiredRxns,1000); % allow only required mets
modelRaven = setParam(modelRaven,'lb',{'bmOUT' 'cordycepinOUT'},0);
modelRaven = setParam(modelRaven,'ub',{'bmOUT' },1000);
modelRaven = setParam(modelRaven,'obj',{'bmOUT'},1);


%{'gabaIN';'glcIN';'glcnIN'}


model_glucose_ref = setParam(modelRaven,'ub',{'glcIN' 'nh3IN'},2); % allow only required mets, glucose and ammonia
model_glucosamine = setParam(modelRaven,'ub','glcnIN',2); % allow only required mets and glucosamine
model_gaba = setParam(modelRaven,'ub','gabaIN',2); % allow only required mets and gaba


sol_model_glucose_ref = solveLP(model_glucose_ref,1);
sol_model_glucosamine = solveLP(model_glucosamine,1);
sol_model_gaba = solveLP(model_gaba,1);


followChanged(modelRaven,sol_model_glucosamine.x,sol_model_glucose_ref.x, 100);
followChanged(modelRaven,sol_model_gaba.x,sol_model_glucose_ref.x, 100);

uniqueDiff = {'r0555';'r0088';'r0031';'r0618';'r0619';'r0620';'r0308';'r0133';'r0158';'r0242';'r0251';'r0153';'r0348';'r0349';'r0399';'r0129';'r0152';'r0342';'r0343';'r0344';'r0345';'r0346';'3amp';'cordycepin';'r0387';'r0557';'r0391';'r0435';'r0553';'r0004';'r0373';'r0541';'r0006';'r0340';'r0341';'r0350';'r1217';'r1226';'r1230';'r1231';'r1316';'r1317';'r0150';'r0347';'r0336';'r0030';'cordycepinOUT';'glcIN';'nh3IN';'nr0010';'r0374';'r1164';'r1189';'r1190';'r1191';'r1238';'r1240'};
rmRxns = setdiff(modelRaven.rxns,uniqueDiff);
subNetwork = removeReactions(modelRaven,rmRxns,true,true);
subNetwork.rxnConfidenceScores = constructEquations(subNetwork);

exportToExcelFormat(subNetwork,'diffsubNetwork.xlsx');


%To identify active and sometime active reactions

model_fixed_glucose = setParam(modelRaven,'ub',{'glcIN'},2); % fixed glucose 
model_fixed_glucose = setParam(model_fixed_glucose,'ub',{'nh3IN'},1000); % allow ammonia

cobraModel_glucose_ref = ravenCobraWrapper(model_glucose_ref);

cobraModel_glucosamine = ravenCobraWrapper(model_glucosamine);
cobraModel_gaba = ravenCobraWrapper(model_gaba);

[grRatio, grRateKO_nh3, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(cobraModel_glucose_ref, 'FBA');
[grRatio, grRateKO_glucosamine, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(cobraModel_glucosamine, 'FBA');
[grRatio, grRateKO_gaba, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(cobraModel_gaba, 'FBA');




model_fixed_glucose = setParam(modelRaven,'ub',{'glcIN'},1000); % fixed glucose 
model_fixed_glucose = setParam(model_fixed_glucose,'ub',{'nh3IN'},1000); % allow ammonia
model_fixed_glucose = setParam(model_fixed_glucose,'ub',{'bmOUT' 'cordycepinOUT'},1000);
model_fixed_glucose = setParam(model_fixed_glucose,'obj',{'bmOUT'},1);

sol = 

cobraModel = readCbModel('C:\Users\Dell\Documents\GitHub\cordyceps_model\ModelFiles\xml\modelBigger.xml');
cobraModel = changeB (cobraModel,'ub',{'glcIN'},1000); % fixed glucose 
model_fixed_glucose = setParam(model_fixed_glucose,'ub',{'nh3IN'},1000); % allow ammonia
model_fixed_glucose = setParam(model_fixed_glucose,'ub',{'bmOUT' 'cordycepinOUT'},1000);
model_fixed_glucose = setParam(model_fixed_glucose,'obj',{'bmOUT'},1);
cobraModel = changeObjective(cobraModel,'bmOUT',1);
sol = optimizeCbModel(cobraModel)
solveCobraLP
[growthRates,shadowPrices1,shadowPrices2] = phenotypePhasePlane(cobraModel,'nh3IN','glcIN', 50, 10, 10)
[growthRates] = phenotypePhasePlane(cobraModel, 'nh3IN', 'glcIN', 5, 20, 20)


