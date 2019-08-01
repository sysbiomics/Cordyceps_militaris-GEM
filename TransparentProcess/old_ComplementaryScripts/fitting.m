%% (Transparent) Parameter Fitting process
% Purposed: This script developed for fitting parameter to Cordyceps militaris model.
% 
% Written by Nachon Raethong, 14-DEC-2018
% Re-upload, 31-JAN-2019
% 
%% WORKSPACE
cd 'C:\Users\Nachonase\Documents\GitHub\cordyceps_model\TransparentProcess';

%% constraint for aerobic growth in minimal media
model = importModel('ModelFiles\xml\model.xml');
model=setParam(model,'lb','cmt_ex_c_glucose_c_e',-0.15);
model=setParam(model,'lb',{'cmt_ex_cordycepin_c_e'},0);
model=setParam(model,'obj','cmt_ex_growth_c',1);
sol = solveLP(model);
printFluxes(model,sol.x,true);

carbons = {'cmt_ex_c_glucose_c_e','cmt_ex_c_fructose_c_e','cmt_ex_c_arabinose_c_e',...
    'cmt_ex_c_xylose_c_e','cmt_ex_c_sucrose_c_e'};
model=setParam(model,'lb',carbons,0);

cordycepin = [0.00234,0.00114,0.00029,0.00060,0.00201];
uptake = [-0.15,-0.125,-0.1,-0.05,-0.08];
for i = 1:numel(uptake)
    model=setParam(model,'lb',carbons,0);
    model=setParam(model,'lb',carbons(i),uptake(i));
    model=setParam(model,'lb',{'cmt_ex_growth_c'},0);
    sugar = carbons(i);
    sugar = replace(sugar,'cmt_ex_c_','');
    sugar = replace(sugar,'_c_e','');
    model=setParam(model,'obj',{'cmt_ex_growth_c'},1);
    sol = solveLP(model);
    fprintf(['Yield of biomass (umax) on ' sugar{1} ' is ' num2str(sol.f*-1) ' /h\n']);
    %model=setParam(model,'lb',{'cmt_ex_cordycepin_c_e'},cordycepin(i));
    %solc = solveLP(model);
    %model=setParam(model,'obj',{'cmt_ex_cordycepin_c_e'},1);
    %solc=solveLP(model);
    %fprintf(['Yield of biomass with cordycepin is ' num2str((solc.f*-1)) ' /h\n']);
   % fprintf(['Yield of cordycepin on biomass is ' num2str(cordycepin(i)/(solc.f*-1)) ' mol/gDW\n']);
end


