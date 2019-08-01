function requiredSubstrates = requiredSubstrates(model)
uptake.Rxns = model.rxns(strncmp(model.rxnNames,'uptake ',7));
uptake.RxnNames = model.rxnNames(strncmp(model.rxnNames,'uptake ',7));
model = setParam(model,'ub',uptake.Rxns,1000);
model = setParam(model,'obj',{'bmOUT'},1);
uptake.umax = cell(numel(uptake.Rxns),1);
uptake.flux = cell(numel(uptake.Rxns),1);
for i = 1:numel(uptake.Rxns)
    modelS = setParam(model,'eq',uptake.Rxns{i},0);
    modelS = setParam(modelS,'obj','bmOUT',1);
    solS = solveLP(modelS,1);
    uptake.umax{i} = solS.f;
    uptake.flux{i} = solS.x;
end


requiredSubstrates = uptake;
end

