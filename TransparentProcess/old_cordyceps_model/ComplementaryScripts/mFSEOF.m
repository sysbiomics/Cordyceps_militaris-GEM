% This is the script used to perform FSEOF 
% Cordycepin production is the target. 
% Input model is iNR818

% Nachon Raethong, 2019-02-24


clear all;
cd 'C:\Users\Nachonase\Documents\GitHub\cordyceps_model';
load 'ModelFiles\mat\model.mat';

% set empty subSystems to be blank 
% usig Hao's script

for i = 1:numel(model.subSystems)
    if ~iscellstr(model.subSystems{i, 1})
        model.subSystems{i, 1} = {'blank'};
    else
    end
end

% constraint the model by allowed only metabolites presence in culture medium
uptakeRxns = model.rxns(strncmp(model.rxnNames,'uptake ',7));
uptakeRxnNames = model.rxnNames(strncmp(model.rxnNames,'uptake ',7));
requiredRxns = {'o2IN' 'piIN' 'slfIN' 'nh3IN'};
model = setParam(model,'ub',uptakeRxns,0);
model = setParam(model,'ub',requiredRxns,1000);
model = setParam(model,'lb',{'bmOUT' 'cordycepinOUT'},0);
model = setParam(model,'ub',{'bmOUT' 'cordycepinOUT'},1000);
% growth (biomass[e]) is an objective.
model = setParam(model,'obj',{'bmOUT'},1);

% FSEOF simulation was performed under three different carbon sources including glucose, xylose and sucrose.
model_gluIN = setParam(model,'ub','glcIN',1);
model_xylIN = setParam(model,'ub','xylIN',1);
model_sucIN = setParam(model,'ub','sucIN',1);

% cordycepin was set to be a target.
targeRxns_gluIN = FSEOF(model_gluIN,'bmOUT','cordycepinOUT');
targeRxns_xylIN = FSEOF(model_xylIN,'bmOUT','cordycepinOUT');
targeRxns_sucIN = FSEOF(model_sucIN,'bmOUT','cordycepinOUT');

save('C:\Users\Nachonase\Documents\GitHub\cordyceps_model\FSEOF\targeRxns_gluIN.mat','targeRxns_gluIN');
save('C:\Users\Nachonase\Documents\GitHub\cordyceps_model\FSEOF\targeRxns_xylIN.mat','targeRxns_xylIN');
save('C:\Users\Nachonase\Documents\GitHub\cordyceps_model\FSEOF\targeRxns_sucIN.mat','targeRxns_sucIN');







