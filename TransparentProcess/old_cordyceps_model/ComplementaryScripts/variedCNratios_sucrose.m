sugars = {'glcIN' 'fruIN' 'xylIN' 'arabIN' 'sucIN'};
model = setParam(modelRaven,'ub',sugars,0);
model = setParam(model,'ub','cordycepinOUT',1000);
model = setParam(model,'lb','cordycepinOUT',0);
model = setParam(model,'obj','bmOUT',1);

variedN = [0.05:0.04:1.05];
CNratios = cell(numel(variedN),1);
umax = cell(numel(variedN),1);
cordycepin = cell(numel(variedN),1);
sucrose = cell(numel(variedN),1);
ammonia = cell(numel(variedN),1);
%glucose
for i = 1:numel(variedN)
    variedC = 1.1 - variedN(i);
    model = setParam(model,'ub','sucIN',variedC);
    model = setParam(model,'ub','nh3IN',variedN(i));
    sol = solveLP(model,1);
    CNratios{i} = num2str((variedC*12)/(variedN(i)*1));
    umax{i} = num2str((sol.f*-1)/1.1);
    model2 = setParam(model,'lb','bmOUT',((sol.f*-1)/1.1));
    model2 = setParam(model2,'obj','cordycepinOUT',1);
    sol2 = solveLP(model2,1);
    cordycepin{i} = num2str(sol2.f*-1);
    sucrose{i} = variedC;
    ammonia{i} = variedN(i);

end

growthProfile_sucIN = table(CNratios,umax,cordycepin,sucrose,ammonia);
