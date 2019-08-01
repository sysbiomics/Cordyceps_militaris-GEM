function newModel = removeUnusedComp(model,UnusedComp)
% This is the function used to remove unused compartment 
%
% INPUT:
%		model = model
%		UnusedComp = compartment id of which compartment do you want to remove
% OUTPUT:
%		newModel = the model without that compartment
%
% Usage:
% 			newModel = removeUnusedComp(model,UnusedComp)
%
% Written by Nachon Raethong, 2019-02-24
% Updated by Nachon Raethong, 2019-08-02
%
%

model = copyToComps(model,UnusedComp,model.rxns(1));
newModel = removeReactions(model,model.rxns(numel(model.rxns)),true,true,true);
end

