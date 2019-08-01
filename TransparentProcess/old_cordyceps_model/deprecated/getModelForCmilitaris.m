%% Reconstruction process
% Purposed: Here is a complete script used for reconstructing a metabolic model of Cordyceps militaris. The % reconstruction process was carried through the following three steps involving data preparation, draft 
% reconstruction and model curation. With its availability of whole-genome sequences, top-down 
% bioinformatics (TDB) approach was mainly applied during the first and second steps to reconstruct a draft 
% metabolic network for C. militaris using RAVEN 2.0. Afterward bottom-up biochemical (BUB) approach 
% was implemented to the last step for enhancing network connectivity, feasibility and capability according % to experimental supports and further gaining a better coverage of the metabolism for C. militaris. This is an 
% iterative process which ends up with a curated metabolic network that can represent the cellular 
% physiology of C. militaris such as growth and cordycepin productivity.
% 
% Written by Nachon Raethong, 31-JAN-2019
% 
%% WORKSPACE
cd 'C:\Users\Dell\Documents\GitHub\cordyceps_model';
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
[~, newMets]=xlsread('ComplementaryData\supplementary_final.xlsx','Table S1');
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
 
 
[~, textData1]=xlsread('ComplementaryData\supplementary_final.xlsx','Table S2');
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
  
[~, textData2]=xlsread('ComplementaryData\supplementary_final.xlsx','Table S3.1');
rxnToAdd = struct();
rxnToAdd.rxns = textData2(4:end,2);
rxnToAdd.grRules = textData2(4:end,10);
curatedModel_v0 = addRxnsGenesMets(cmtDraftFromPenicillium,reduced,rxnToAdd.rxns,...
    rxnToAdd.grRules,'additional rxns from Penicillium model with gene annotation from alternative templates',2);
 

%% STEP 3: MODEL CURATION
% Addition of biomass reactions and some new reactions 
%	New reactions which excluded from iAL1006 were manually formulated and added with their 
% annotated genes, associated reactions and assigned compartments by the following steps. 
% First, new formulated reactions listed in Table S4 were manually added by addRxns function to the 
% template model for check mass balance of each reaction. 
 
[~, SheetS]=xlsread('ComplementaryData\supplementary_final.xlsx','Table S4');
cmtBiomass = struct();
cmtBiomass.rxns = SheetS(2:end,1);
cmtBiomass.rxnNames = SheetS(2:end,2);
cmtBiomass.equations = SheetS(2:end,3);
optimizedBiomass = addRxns(reduced,cmtBiomass,2,'c',true);
optimizedBiomass = setParam(optimizedBiomass,'lb',cmtBiomass.rxns,0);
optimizedBiomass = setParam(optimizedBiomass,'ub',cmtBiomass.rxns,1000);
optimizedBiomass = setParam(optimizedBiomass,'eq','r1463',0);
optimizedBiomass = setParam(optimizedBiomass,'eq','r1466',0);
optimizedBiomass = setParam(optimizedBiomass,'lb','cordycepinOUT',0);
optimizedBiomass = setParam(optimizedBiomass,'eq','r1467',0);
optimizedBiomass = setParam(optimizedBiomass,'eq','nr0007',1);
optimizedBiomass = setParam(optimizedBiomass,'obj','bmOUT',1);
 
sol = solveLP(optimizedBiomass);
printFluxes(optimizedBiomass, sol.x);

% Then, biomass and new reactions listed in Table S4 were added to the curated metabolic network of C. militaris 
% by recruited from the template model that included these reactions (namely, optimizedBiomass).

cmtBiomass.grRules = SheetS(2:end,4);
curatedModel_v1 = addRxnsGenesMets(curatedModel_v0 ,optimizedBiomass,cmtBiomass.rxns,...
    cmtBiomass.grRules,'additional rxns for formulating biomass',2);

% To allow flux distribution from substrates such as glucose flow to demand
% objective functions such as growth. Some exchange reactions and required
% sponteous reactions as listed in Table S5 were introduced to the final model.

[~, Sheet8]=xlsread('ComplementaryData\supplementary_final.xlsx','Table S5');
nonGENErxns = Sheet8(1:end,1);

curatedModel_v2 = addRxnsGenesMets(curatedModel_v1,optimizedBiomass,nonGENErxns,false,...
    'nonGENErxns.rxns for modeling',1);
 
% Feasibility, functionality, capability and validity of the final
% model were checked by maximizing growth and constrainted lower bound of cordycepin production 
% to be 0.0001 to allow the production of cordycepin during growth on different substrates.  
curatedModel = curatedModel_v2;
curatedModel = setParam(curatedModel,'lb','cordycepinOUT',0.0001);
curatedModel = setParam(curatedModel,'eq','r0181',0);
curatedModel = setParam(curatedModel,'eq','r0019',0);
curatedModel = setParam(curatedModel,'obj',{'bmOUT'},1);
sol = solveLP(curatedModel);
sugars = {'glcIN' 'fruIN' 'xylIN' 'arabIN' 'sucIN'};
uptake = [0.15,0.125,0.054,0.11,0.07499];
sumax(1) = sol;
for i = 1:numel(uptake)
    model=setParam(curatedModel,'ub',sugars,0);
    model=setParam(model,'ub',sugars(i),uptake(i));
    sugar = sugars(i);
    sol = solveLP(model,1);
    sumax(i) = sol;
    umax = num2str(sol.f*-1);
    fprintf(['Yield of biomass (umax) on ' sugar{1} ' is ' num2str(sol.f*-1) ' /h\n']);
end

% By the end the predicted values are fitted with the particular
% experimental values
%% Finalize model for export to Github
 
curatedModel.id = 'iNR818';
curatedModel.description = 'Genome-scale metabolic model of C militaris by Nachon [31-JAN-2019]';
curatedModel.annotation = [];
exportForGit(curatedModel,'model','C:\Users\Nachonase\Documents\GitHub\cordyceps_model\');
