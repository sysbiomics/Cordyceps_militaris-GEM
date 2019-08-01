%% Reconstruction process
% Purposed: Here is a complete script used for reconstructing a metabolic model of Cordyceps militaris.  
% The reconstruction process was carried through the following 4 steps. 
% This is an iterative process which ends up with a validated GEM 
% that can represent C. militaris growth on different carbon sources.
% 
% Written by Nachon Raethong, 31-JAN-2019
% Updated by Nachon Raethong, 24-JULY-2019

%% WORKSPACE
cd 'C:\Users\Dell\Documents\GitHub\Cordyceps_militaris-GEM';
%% STEP 1: DATA PREPARATION
% 1.1 Data collection
% 	Whole-proteome sequences of C. militaris and template species as well as the model files for 
% Penicillium chrysogenum Wisconsin 54-1255 (iAL1006), Aspergillus nidulans FGSC A4 (iHD666) and
% Neurospora crassa (iJDZ836) were collected and keep in /ComplementaryData/ folder. Notably, gene
% identifiers of these proteome sequences must be same with gene identifiers used in template models.

peniModel = importModel('ComplementaryData\iAL1006 v1.00.xml');
anidModel = importModel('ComplementaryData\ANIDULANS.xml');
neurModel = importModel('ComplementaryData\iJDZ808.xml');


% 1.2 Addition of new metabolites to the templates network  
% 	New metabolites listed in Table S1 were introduced to the templates network by addMets function.
[~, newMets]=xlsread('ComplementaryData\supplementary.xlsx','Table S1');
metsToAdd = struct();
metsToAdd.mets = newMets(2:end,1);
metsToAdd.metNames = newMets(2:end,2);
metsToAdd.metFormulas = newMets(2:end,3);
metsToAdd.compartments = 'c';
peniModel_newMets=addMets(peniModel,metsToAdd);

% 1.3 Metabolite naming
% 	The metabolite names used by different template models may introduce a lot of duplicated 
% metabolites and reactions during the reconstruction process. So, metabolites names used in these template
% models were reorganized by mapping metabolite identifiers obtained from different databases (e.g. 
% MetaCyc, KEGG or CHEBI) to the same metabolite names thereby resulted in compatible metabolite 
% namespaces. Overall updated metabolite names and their identifiers are listed in Table S2.
 
 
[~, textData1]=xlsread('ComplementaryData\supplementary.xlsx','Table S2');
metNames = struct();
metNames.old = textData1(4:end,2);
metNames.new = textData1(4:end,3);
[a, b]=ismember(peniModel_newMets.metNames,metNames.old);
I=find(a);
peniModel_newMets.metNames(I)=metNames.new(b(I));
reduced=contractModel(peniModel_newMets);
reduced = setParam(reduced,'ub',{'glcIN' 'etohIN'},[1 0]);
reduced = setParam(reduced,'obj',{'bmOUT'},1);
sol = solveLP(reduced);
printFluxes(reduced, sol.x, true);
idex_cordycepinPermease = find(ismember(reduced.rxns,'r1468'));
reduced.rxns(idex_cordycepinPermease) = {'nr0010'};
reduced.rxnNames(idex_cordycepinPermease) = {'cordycepin permease'};
idex_cordycepinProduction = find(ismember(reduced.rxns,'penartOUT'));
reduced.rxns(idex_cordycepinProduction) = {'cordycepinOUT'};
reduced.rxnNames(idex_cordycepinProduction) = {'production of cordycepin'};
 

%% STEP 2: DRAFT RECONSTRUCTION
% 2.1 Identification of orthologous proteins 
%	A set of orthologous proteins between C. militaris and template model proteins was obtained by 
% sequence-alignment analysis using BLASTp.
load 'ComplementaryData\blastedCordycepsAnidulans.mat';
load 'ComplementaryData\blastedCordycepsPenicillium.mat';
load 'ComplementaryData\blastedCordycepsNcrassa.mat';
 
% 2.2 Get draft networks from protein orthology inference
% 	 Orthologous proteins were mapped to template models using getModelFromHomology function. 
cmtDraftFromAnidulans = getModelFromHomology({anidModel},blastedCordycepsAnidulans,'Cordyceps',{},2,false,10^-5,100);
cmtDraftFromAnidulans.id = 'cmtDraftFromAnidulans';
cmtDraftFromNcrassa = getModelFromHomology({neurModel},blastedCordycepsNcrassa,'Cordyceps',{},2,false,10^-5,100);
cmtDraftFromNcrassa.id = 'cmtDraftFromNcrassa';
cmtDraftFromPenicillium = getModelFromHomology({reduced},blastedCordycepsPenicillium,'Cordyceps',{},2,false,10^-5,100);
cmtDraftFromPenicillium.id = 'cmtDraftFromPenicillium';
 

% 2.3 Enhanced metabolic coverage of the main constructed network by alternative networks
% 	In this case, iAL1006 is the first model reconstructed by RAVEN and intensively used for RAVEN 2.0 
% tutorial, so the constructed draft initiated from Penicillium model (iAL1006) is promisingly used as a main 
% network for C. militaris. In particular other drafts created from iHD666 and iJDZ836 as well as an earlier 
% genome-scale metabolic network of C. militaris (iWV1170) were used as alternative networks supported 
% for enhancing metabolic coverage of the main network. For instance, by checking the EC numbers, some 
% new reactions were added from iAL1006 with annotated genes from alternative networks as listed in Table 
% S3.1 using addGenesMetsRxns function. The comparison is shown in Table S3.2.
  
[~, textData2]=xlsread('ComplementaryData\supplementary.xlsx','Table S3.1');
rxnToAdd = struct();
rxnToAdd.rxns = textData2(4:end,2);
rxnToAdd.grRules = textData2(4:end,10);
curatedModel_v0 = addRxnsGenesMets(cmtDraftFromPenicillium,reduced,rxnToAdd.rxns,...
    rxnToAdd.grRules,'additional rxns from Penicillium model with gene annotation from alternative templates',2);
 

%% STEP 3: MODEL CURATION
% Addition of biomass reactions and some new reactions 
% New reactions were manually formulated and added with their 
% annotated genes, associated reactions compartments and assigned subSystems. 
% First, new formulated reactions listed in Table S4 were manually added by addRxns function to the 
% template model for check mass balance of each reaction. 

[~, SheetS]=xlsread('ComplementaryData\supplementary.xlsx','Table S4');
cmtBiomass = struct();
cmtBiomass.rxns = SheetS(2:end,1);
cmtBiomass.rxnNames = SheetS(2:end,2);
cmtBiomass.equations = SheetS(2:end,3);
cmtBiomass.grRules = SheetS(2:end,4);

% try to add C. militaris biomass rxns to reduced model of iAL1006 for
% checking validity of these constructed rxns before add to C. militaris
% model
optimizedBiomass = addRxns(reduced,cmtBiomass,2,'c',true,true);
optimizedBiomass = setParam(optimizedBiomass,'lb',cmtBiomass.rxns,0);
optimizedBiomass = setParam(optimizedBiomass,'ub',cmtBiomass.rxns,1000);
optimizedBiomass = setParam(optimizedBiomass,'eq','r1463',0); %block penicillium biomass reaction {'0.04 amino acid pool[c] + 104 ATP[c] + 0.25 cellwall[c] + 0.0001 cofactors[c] + 0.01 deoxyribonucleic acids[c] + 0.005 glycerides and free fatty acids[c] + 104 H2O[c] + 1e-05 other lipids[c] + 0.45 protein[c] + 0.08 ribonucleic acids[c] + 0.01 sterol esters[c] + 0.035 phospholipids[e] => 104 ADP[c] + biomass[c] + 104 phosphate[c]'}
optimizedBiomass = setParam(optimizedBiomass,'eq','r1466',0); %block penicillium ATP reaction {'4.14 ATP[c] + 4.14 H2O[c] => 4.14 ADP[c] + 4.14 phosphate[c]'}
optimizedBiomass = setParam(optimizedBiomass,'eq','r1467',0); %block penicillium penicilin reaction {'ATP[p] + H2O[p] + isopenicillin N[p] => AMP[p] + artificial penicillin[p] + diphosphate[p] + L-2-aminoadipate[p]'}
optimizedBiomass = setParam(optimizedBiomass,'lb','cordycepinOUT',0);
optimizedBiomass = setParam(optimizedBiomass,'ub',{'proteinIN','dnaIN','rnaIN','lipidIN','carbohydrateIN'},[0 0 0 0 0]);
optimizedBiomass = setParam(optimizedBiomass,'obj','bmOUT',1);
 
sol = solveLP(optimizedBiomass);
printFluxes(optimizedBiomass, sol.x);

% If reduced model of iAL1006 can grow by C. militaris biomass rxns, we
% will add these constructed biomass rxns to the curated metabolic network of C. militaris 
% by recruited from the reduced model of iAL1006 that included these
% reactions (namely, optimizedBiomass).
curatedModel_v1 = addRxnsGenesMets(curatedModel_v0 ,optimizedBiomass,cmtBiomass.rxns,...
    cmtBiomass.grRules,'additional rxns for formulating biomass',2);

% To allow flux distribution from substrates such as glucose flow to demand
% objective functions such as growth. Some exchange reactions and required
% sponteous reactions as listed in Table S5 were introduced to the final model.

[~, Sheet8]=xlsread('ComplementaryData\supplementary.xlsx','Table S5');
nonGENErxns = Sheet8(1:end,1);

curatedModel_v2 = addRxnsGenesMets(curatedModel_v1,optimizedBiomass,nonGENErxns,false,...
    'nonGENErxns.rxns for modeling',1);
 
% Feasibility, functionality, capability and validity of the final
% model were checked by maximizing growth and constrainted lower bound of cordycepin production 
curatedModel = removeUnusedComp(curatedModel_v2,'b');
curatedModel = setParam(curatedModel,'obj',{'bmOUT'},1);

%To expand curatedModel with iWV1346 
cd 'C:\Users\Dell\Documents\GitHub\Cordyceps_militaris-GEM';
aor = importModel('ComplementaryData\iWV1346.xml');
ccm = curatedModel;
%blastedCordycepsAoryzae = getBlast('Cordyceps','ComplementaryData\Cmilitaris_protein.fasta','iWV1346','ComplementaryData\A_oryzae.faa');
%save('ComplementaryData\blastedCordycepsAoryzae.mat','blastedCordycepsAoryzae');
aorModel = aor;
[~, newMets]=xlsread('ComplementaryData\supplementary.xlsx','Table S1.2');
metsToAdd = struct();
metsToAdd.mets = newMets(2:end,1);
metsToAdd.metNames = newMets(2:end,2);
metsToAdd.metFormulas = newMets(2:end,3);
metsToAdd.compartments = 'c';
aorModel_newMets=addMets(aorModel,metsToAdd);


[~, textData1]=xlsread('ComplementaryData\supplementary.xlsx','Table S2');
metNames = struct();
metNames.old = textData1(4:end,2);
metNames.new = textData1(4:end,3);
[a, b]=ismember(aorModel_newMets.metNames,metNames.old);
I=find(a);
aorModel_newMets.metNames(I)=metNames.new(b(I));
aorModel_newMets.metNames = lower(aorModel_newMets.metNames);
reduced=contractModel(aorModel_newMets);
load 'ComplementaryData\blastedCordycepsAoryzae.mat';
cmtDraftFromAoryzae = getModelFromHomology({reduced},blastedCordycepsAoryzae,'Cordyceps',{},2,false,10^-5,100);
for i = 1:numel(cmtDraftFromAoryzae.rxns)
    cmtDraftFromAoryzae.lb(i) = 0;
    cmtDraftFromAoryzae.ub(i) = 0;
end

merged = addRxnsGenesMets(ccm,cmtDraftFromAoryzae,cmtDraftFromAoryzae.rxns,true);
new = contractModel(merged);


%To expand the current C. militaris model (namely, new) with KEGG models

load 'ComplementaryData\cmtKEGG.mat';
pre_cmtKEGG = removeMets(cmtKEGG,'H+',true,true,false);
pre_cmtKEGG.metNames = lower(pre_cmtKEGG.metNames);
[pre_cmtKEGG, removedRxns] = removeBadRxns(pre_cmtKEGG);
pre_cmtKEGG = removeMets(pre_cmtKEGG,'udp (g10619)',true,true,false);
pre_cmtKEGG = removeMets(pre_cmtKEGG,'e-',true,true,false);

[~, textData1]=xlsread('ComplementaryData\supplementary.xlsx','Table S2');
metNames = struct();
metNames.old = textData1(4:end,2);
metNames.new = textData1(4:end,3);
[a, b]=ismember(pre_cmtKEGG.metNames,metNames.old);
I=find(a);
pre_cmtKEGG.metNames(I)=metNames.new(b(I));
pre_cmtKEGG = copyToComps(pre_cmtKEGG,'c',pre_cmtKEGG.rxns,true,'Cytosol');
reducedKEGG=contractModel(pre_cmtKEGG);
for i = 1:numel(reducedKEGG.rxns)
    reducedKEGG.rev(i) = 0;
    reducedKEGG.lb(i) = 0;
    reducedKEGG.ub(i) = 0;
end

merge_ccm_aor_kegg = addRxnsGenesMets(new,reducedKEGG,reducedKEGG.rxns,true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%remove bad rxns
[~, badRxns]=xlsread('ComplementaryData\supplementary.xlsx','Table S6');
reactionsToRemove = badRxns(1:end,1);
model = removeReactions(merge_ccm_aor_kegg,reactionsToRemove,true,true);
[~, textData9]=xlsread('ComplementaryData\supplementary.xlsx','Table S7');
subSystem = struct();
subSystem.rxns = textData9(2:end,1);
subSystem.new = textData9(2:end,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newSubSystemModel = contractModel(model);
[a, b]=ismember(newSubSystemModel.rxns,subSystem.rxns);
I=find(a);
newSubSystemModel.rxnNotes(I)=subSystem.new(b(I));
modelBigger = contractModel(newSubSystemModel);
[~, textData10]=xlsread('ComplementaryData\supplementary.xlsx','Table S8');
modelBigger.rxnReferences = textData10(2:end,2);

%% STEP 4: MODEL VALIDATION
% The model validation was performed against experimentation.
validateModel = modelBigger;
validateModel = setParam(validateModel,'lb',{'bmOUT','cordycepinOUT'},0);
validateModel = setParam(validateModel,'eq',{'matp'},1);
validateModel = setParam(validateModel,'ub',{'bmOUT','cordycepinOUT'},1000);
validateModel = setParam(validateModel,'obj',{'bmOUT'},1);

% The following information was reported in Table 2
sugars = {'glcIN' 'fruIN' 'arabIN' 'xylIN' 'sucIN'};
uptake = [0.1448,0.1472,0.1074,0.0681,0.0815]; %observed from experiment
sumax(1) = sol;

for i = 1:numel(uptake)
    model=setParam(validateModel,'ub',sugars,0);
    model=setParam(model,'ub',sugars(i),uptake(i));
    sugar = sugars(i);
    sol = solveLP(model,1);
    sumax(i) = sol;
    umax = num2str(sol.f*-1);
    fprintf([num2str(sol.f*-1) '\n']);
end

% If the model prediction matched the experimentation with the Error rates below
% 5%, then we will export the model for the latter analysis.

validateModel.id = 'iNR1320';
validateModel.description = 'C. militaris GEM by Nachon Raethong [24-JULY-2019]';
validateModel.annotation = [];
exportForGit(validateModel,'model','C:\Users\Dell\Documents\GitHub\Cordyceps_militaris-GEM\');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%