cd 'C:\Users\Dell\Documents\GitHub\cordyceps_model';
load('C:\Users\Dell\Documents\GitHub\cordyceps_model\ModelFiles\mat\modelBigger.mat');
modelBigger = model;
sol = solveLP(modelBigger,1);
fprintf(['Yield of biomass (umax) is '  num2str(sol.f*-1) ' /h\n']);


modelRaven = modelBigger;
%uptakeRxns = {'glcIN' 'nh3IN' 'sucIN' 'glnIN' 'chitIN'};
uptakeRxns = modelRaven.rxns(strncmp(modelRaven.rxnNames,'uptake ',7));
uptakeRxnNames = modelRaven.rxnNames(strncmp(modelRaven.rxnNames,'uptake ',7));
requiredRxns = {'o2IN' 'piIN' 'slfIN'};

modelRaven = setParam(modelRaven,'ub',uptakeRxns,0); % block all uptake
modelRaven = setParam(modelRaven,'ub',requiredRxns,1000); % allow only required mets
modelRaven = setParam(modelRaven,'lb',{'bmOUT' 'cordycepinOUT'},0);
modelRaven = setParam(modelRaven,'ub',{'bmOUT' 'cordycepinOUT'},1000);
modelRaven = setParam(modelRaven,'obj',{'bmOUT'},1);

mVersatile_I = cell(numel(uptakeRxns),1);
solVersatile_I= cell(numel(uptakeRxns),1);
model_mVersatile_I = modelRaven; % allow only required mets 
for i = 1:numel(uptakeRxns)
    model = model_mVersatile_I;
    model_i = setParam(model,'ub',uptakeRxns{i},2);
    sol = solveLP(model_i,1);
    mVersatile_I{i} = (sol.f*-1);
    solVersatile_I{i} = sol;
end

solVersatile_II= cell(numel(uptakeRxns),1);
mVersatile_II = cell(numel(uptakeRxns),1);
model_mVersatile_II = setParam(modelRaven,'ub','glcIN',1); % allow only required mets and glcIN
for i = 1:numel(uptakeRxns)
    model = model_mVersatile_II;
    model_i = setParam(model,'ub',uptakeRxns{i},1);
    sol = solveLP(model_i,1);
    mVersatile_II{i} = (sol.f*-1);
    solVersatile_II{i} = sol;
end


solVersatile_III = cell(numel(uptakeRxns),1);
mVersatile_III = cell(numel(uptakeRxns),1);
model_mVersatile_III = setParam(modelRaven,'ub','nh3IN',1); % allow only required mets and nh4IN
for i = 1:numel(uptakeRxns)
    model = model_mVersatile_III;
    model_i = setParam(model,'ub',uptakeRxns{i},1);
    sol = solveLP(model_i,1);
    mVersatile_III{i} = (sol.f*-1);
    solVersatile_III{i} = sol;
end


variedResult = table(uptakeRxns,uptakeRxnNames,mVersatile_I,mVersatile_II,mVersatile_III);