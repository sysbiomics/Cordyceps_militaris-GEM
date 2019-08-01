cd 'C:\Users\Dell\Documents\GitHub\cordyceps_model';
load('C:\Users\Dell\Documents\GitHub\cordyceps_model\ModelFiles\mat\modelBigger.mat');
validateModel = curatedModel;
validateModel = setParam(validateModel,'lb',{'bmOUT','cordycepinOUT'},0);
validateModel = setParam(validateModel,'eq',{'matp'},1);
validateModel = setParam(validateModel,'ub',{'bmOUT','cordycepinOUT'},1000);
validateModel = setParam(validateModel,'obj',{'bmOUT'},1);

sol = solveLP(validateModel,1);
printFluxes(validateModel,sol.x,true);
fprintf(['Yield of biomass (umax) is '  num2str(sol.f*-1) ' /h\n']);

sugars = {'glcIN' 'fruIN' 'arabIN' 'xylIN' 'sucIN'};
uptake = [0.15499,0.15,0.1,0.07499,0.09];
sumax(1) = sol;



for i = 1:numel(uptake)
    model=setParam(validateModel,'ub',sugars,0);
    model=setParam(model,'ub',sugars(i),uptake(i));
    sugar = sugars(i);
    sol = solveLP(model,1);
    sumax(i) = sol;
    umax = num2str(sol.f*-1);
    %fprintf(['Maximum yield of biomass (umax) on ' sugar{1} ' is ' num2str(sol.f*-1) ' /h\n']);
    fprintf([num2str(sol.f*-1) '\n']);
end
