%% POPCORNtutorial.m
%% Objective: to use POPCORN
% This script was written for showing how to use POPCORN by the example
% case on the optimization of C:N ratio for fast growth and cordycepin 
% production in Cordyceps militaris using iNR1320.
%
% Written by Nachon Raethong, 08-AUG-2019
%

cd 'C:\Users\Dell\Documents\GitHub\Cordyceps_militaris-GEM';
% load iNR1320
load('C:\Users\Dell\Documents\GitHub\Cordyceps_militaris-GEM\ModelFiles\mat\model.mat');

% prepare the model for POPCORN such as set up the lowwer and upper bounds
% of the biomass and target production reactions
POPCORNmodel = setParam(model,'ub',{'bmOUT';'cordycepinOUT'},[1000 1000]);
POPCORNmodel = setParam(POPCORNmodel,'lb',{'bmOUT';'cordycepinOUT'},[0 0]);

% apply POPCORN by the following arguments
growthProfile = POPCORN(POPCORNmodel,'bmOUT','cordycepinOUT','glcIN','nh3IN',6,1);

