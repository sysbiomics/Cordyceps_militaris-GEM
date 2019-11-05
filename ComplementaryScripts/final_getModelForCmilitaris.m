%% Reconstruction process
% Purposed: Here is a complete script used for reconstructing a metabolic model of Cordyceps militaris.  
% The reconstruction process was carried through the following 4 steps. 
% This is an iterative process which ends up with a validated GEM 
% that can represent C. militaris growth on different carbon sources.
% 
% Written by Nachon Raethong, 31-JAN-2019
% Updated by Nachon Raethong, 24-JULY-2019
% Updated by Nachon Raethong, 29-OCT-2019

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


%% STEP 4.1: MODEL REFINEMENT
% Remove duplicated reactions
refineModel = modelBigger;
%'R09145_c';'R09144_c';
duplicateRxns = {'R09144_c';'r12';'r13';'r120';'r135';'r158';'r185';'r186';'r188';'r191';'r198';'r236';'r238';'r250';'r274';'r334';'r372';'r374';'r377';'r378';'r379';'r458';'r464';'r495';'r517';'r531';'r559';'r609';'r647';'r666';'r685';'r686';'r714';'r717';'r729';'r749';'r754';'r769';'r786';'r787';'r798';'r799';'r800';'r801';'r803';'r804';'r805';'r807';'r811';'r813';'r816';'r819';'r828';'r830';'r831';'r833';'r835';'r836';'r837';'r845';'r857';'r871';'r888';'r896';'r897';'r899';'r900';'r924';'r929';'r939';'r947';'r958';'r964';'r966';'r969';'r971';'r977';'r1156';'r1157';'r1159';'r1163';'r1279';'r1281';'r1395';'r1459';'r1461';'r1464';'r1466';'r1496';'r1498';'r1500';'r1504';'r1507';'r1509';'r1515';'r1533';'r1575';'r1578';'r1581';'r1582';'r1583';'r1676';'r1685';'r1686';'r1688';'r1689';'r1693';'r1728';'r1734';'r1737';'r1739';'r1749';'r1762';'r1768';'r1769';'r1770';'r1779';'r1792';'r1798';'r1817';'r2098';'r2163';'r2179';'R00009_c';'R00028_c';'R00093_c';'R00103_c';'R00112_c';'R00114_c';'R00139_c';'R00156_c';'R00190_c';'R00197_c';'R00200_c';'R00238_c';'R00239_c';'R00245_c';'R00258_c';'R00268_c';'R00289_c';'R00291_c';'R00299_c';'R00310_c';'R00319_c';'R00330_c';'R00342_c';'R00355_c';'R00369_c';'R00372_c';'R00416_c';'R00426_c';'R00428_c';'R00430_c';'R00472_c';'R00479_c';'R00489_c';'R00513_c';'R00516_c';'R00517_c';'R00570_c';'R00588_c';'R00621_c';'R00658_c';'R00678_c';'R00691_c';'R00699_c';'R00702_c';'R00705_c';'R00706_c';'R00707_c';'R00708_c';'R00714_c';'R00715_c';'R00720_c';'R00722_c';'R00736_c';'R00742_c';'R00746_c';'R00755_c';'R00756_c';'R00760_c';'R00762_c';'R00765_c';'R00768_c';'R00774_c';'R00835_c';'R00842_c';'R00848_c';'R00851_c';'R00885_c';'R00891_c';'R00893_c';'R00896_c';'R00908_c';'R00922_c';'R00927_c';'R00935_c';'R00936_c';'R00937_c';'R00939_c';'R00940_c';'R00942_c';'R00945_c';'R00955_c';'R00959_c';'R00962_c';'R00967_c';'R00968_c';'R00970_c';'R00994_c';'R01004_c';'R01015_c';'R01041_c';'R01049_c';'R01057_c';'R01061_c';'R01067_c';'R01071_c';'R01072_c';'R01073_c';'R01082_c';'R01083_c';'R01086_c';'R01098_c';'R01123_c';'R01127_c';'R01135_c';'R01137_c';'R01138_c';'R01163_c';'R01184_c';'R01213_c';'R01224_c';'R01226_c';'R01229_c';'R01230_c';'R01248_c';'R01251_c';'R01252_c';'R01273_c';'R01324_c';'R01325_c';'R01351_c';'R01360_c';'R01364_c';'R01372_c';'R01398_c';'R01404_c';'R01481_c';'R01497_c';'R01518_c';'R01529_c';'R01542_c';'R01548_c';'R01549_c';'R01602_c';'R01641_c';'R01648_c';'R01678_c';'R01698_c';'R01701_c';'R01702_c';'R01752_c';'R01768_c';'R01769_c';'R01773_c';'R01775_c';'R01819_c';'R01830_c';'R01855_c';'R01857_c';'R01858_c';'R01868_c';'R01872_c';'R01880_c';'R01899_c';'R01900_c';'R01914_c';'R01933_c';'R01934_c';'R01936_c';'R01939_c';'R01940_c';'R01954_c';'R01960_c';'R01976_c';'R01978_c';'R01986_c';'R02016_c';'R02017_c';'R02019_c';'R02021_c';'R02024_c';'R02080_c';'R02082_c';'R02085_c';'R02091_c';'R02093_c';'R02094_c';'R02096_c';'R02097_c';'R02098_c';'R02103_c';'R02106_c';'R02107_c';'R02109_c';'R02112_c';'R02133_c';'R02222_c';'R02235_c';'R02236_c';'R02251_c';'R02282_c';'R02283_c';'R02291_c';'R02300_c';'R02301_c';'R02315_c';'R02320_c';'R02326_c';'R02327_c';'R02331_c';'R02332_c';'R02366_c';'R02367_c';'R02371_c';'R02372_c';'R02401_c';'R02410_c';'R02413_c';'R02425_c';'R02433_c';'R02466_c';'R02485_c';'R02487_c';'R02488_c';'R02508_c';'R02519_c';'R02521_c';'R02530_c';'R02536_c';'R02537_c';'R02566_c';'R02570_c';'R02571_c';'R02596_c';'R02619_c';'R02649_c';'R02662_c';'R02665_c';'R02668_c';'R02678_c';'R02695_c';'R02697_c';'R02701_c';'R02703_c';'R02720_c';'R02736_c';'R02740_c';'R02750_c';'R02864_c';'R02872_c';'R02874_c';'R02940_c';'R02957_c';'R02971_c';'R02978_c';'R02984_c';'R03004_c';'R03012_c';'R03024_c';'R03038_c';'R03051_c';'R03139_c';'R03158_c';'R03174_c';'R03180_c';'R03220_c';'R03222_c';'R03236_c';'R03237_c';'R03238_c';'R03239_c';'R03291_c';'R03293_c';'R03300_c';'R03302_c';'R03313_c';'R03316_c';'R03321_c';'R03348_c';'R03355_c';'R03443_c';'R03471_c';'R03530_c';'R03531_c';'R03596_c';'R03635_c';'R03646_c';'R03650_c';'R03652_c';'R03654_c';'R03655_c';'R03657_c';'R03658_c';'R03662_c';'R03663_c';'R03665_c';'R03778_c';'R03815_c';'R03858_c';'R03869_c';'R03909_c';'R03919_c';'R03920_c';'R03921_c';'R03936_c';'R03947_c';'R03968_c';'R03991_c';'R04001_c';'R04007_c';'R04027_c';'R04065_c';'R04097_c';'R04125_c';'R04138_c';'R04144_c';'R04171_c';'R04173_c';'R04178_c';'R04188_c';'R04209_c';'R04225_c';'R04355_c';'R04371_c';'R04378_c';'R04385_c';'R04386_c';'R04391_c';'R04424_c';'R04426_c';'R04428_c';'R04430_c';'R04439_c';'R04440_c';'R04444_c';'R04445_c';'R04506_c';'R04535_c';'R04537_c';'R04544_c';'R04559_c';'R04560_c';'R04568_c';'R04725_c';'R04726_c';'R04742_c';'R04747_c';'R04783_c';'R04862_c';'R04882_c';'R04883_c';'R04888_c';'R04889_c';'R04891_c';'R04892_c';'R04903_c';'R04909_c';'R04928_c';'R04936_c';'R04942_c';'R04944_c';'R04945_c';'R04946_c';'R04952_c';'R04954_c';'R04956_c';'R04957_c';'R04959_c';'R04960_c';'R04962_c';'R04963_c';'R04965_c';'R04967_c';'R04968_c';'R04970_c';'R04972_c';'R04982_c';'R04983_c';'R04988_c';'R04996_c';'R05046_c';'R05050_c';'R05051_c';'R05052_c';'R05064_c';'R05066_c';'R05068_c';'R05069_c';'R05071_c';'R05085_c';'R05231_c';'R05237_c';'R05238_c';'R05286_c';'R05487_c';'R05551_c';'R05576_c';'R05590_c';'R05614_c';'R05639_c';'R05640_c';'R05731_c';'R06004_c';'R06087_c';'R06114_c';'R06154_c';'R06366_c';'R06525_c';'R06526_c';'R06590_c';'R06740_c';'R06941_c';'R07034_c';'R07035_c';'R07104_c';'R07136_c';'R07168_c';'R07363_c';'R07364_c';'R07392_c';'R07411_c';'R07412_c';'R07443_c';'R07460_c';'R07481_c';'R07483_c';'R07494_c';'R07497_c';'R07599_c';'R07600_c';'R07601_c';'R07602_c';'R07603_c';'R07604_c';'R07618_c';'R07888_c';'R07891_c';'R07892_c';'R07895_c';'R07896_c';'R07899_c';'R07934_c';'R07937_c';'R07942_c';'R07950_c';'R07953_c';'R07977_c';'R07978_c';'R07979_c';'R07981_c';'R08090_c';'R08146_c';'R08218_c';'R08221_c';'R08231_c';'R08232_c';'R08235_c';'R08243_c';'R08244_c';'R08282_c';'R08283_c';'R08307_c';'R08363_c';'R08364_c';'R08381_c';'R08639_c';'R08734_c';'R08739_c';'R09076_c';'R09099_c';'R09245_c';'R09246_c';'R09300_c';'R09301_c';'R09372_c';'R09394_c';'R09783_c';'R09951_c';'R09993_c';'R10046_c';'R10052_c';'R10074_c';'R10115_c';'R10119_c';'R10170_c';'R10221_c';'R10619_c';'R10700_c';'R10712_c';'R10907_c';'R10996_c';'R10997_c';'R10998_c';'R11104_c';'R11262_c';'R11316_c';'R11329_c';'R11528_c';'R11529_c';'R11765_c';'R11783_c';'R11893_c';'R11894_c';'R11895_c';'R11906_c';'R02317_c';'R04639_c';'R05048_c';'R06982_c';'R08637_c';'R08698_c';'R11098_c';'R11099_c';'R04241_c';'R00698_c';'R04546_c';'R02918_c';'R00996_c';'R00332_c';'R02090_c';'R05921_c';'R04984_c';'R05641_c';'R10993_c';'R00348_c';'R00269_c';'R08359_c';'R01092_c';'R00006_c';'R00226_c';'R03050_c';'R04672_c';'R04673_c';'R08648_c';'R01068_c';'R01070_c';'R01829_c';'R02568_c';'R00709_c';'R09107_c';'R03661_c';'R03194_c';'R03648_c';'R01699_c';'R03270_c';'R00209_c';'R00014_c';'R05577_c';'R00137_c';'R02058_c';'R01800_c';'R00410_c';'R00248_c';'R00590_c';'R03601_c';'R04859_c';'R11593_c';'R10994_c';'R00830_c';'R00527_c';'R01818_c';'R03444_c';'R02472_c';'R01039_c';'R01799_c';'R05578_c';'R01262_c';'R01687_c';'R03867_c';'R03916_c';'R03970_c';'R03971_c';'R04935_c';'R00508_c';'R00177_c';'R04771_c';'R02164_c';'R02670_c';'R00602_c';'R02517_c';'R07129_c';'R02702_c';'R02909_c';'R03628_c';'R02396_c';'R00259_c';'R04941_c';'R00794_c';'R00796_c';'R07506_c';'R01186_c';'R03629_c';'R04121_c';'R05259_c';'R03098_c';'R04390_c';'R04863_c';'R01209_c';'R04441_c';'R05070_c';'R02111_c';'R03243_c';'R00694_c';'R00734_c';'R00187_c';'R01870_c';'R03660_c';'R00858_c';'R05963_c';'R07809_c';'R07810_c';'R10831_c';'R02199_c';'R10991_c';'R01214_c';'R01090_c';'R03656_c';'R03425_c';'R03033_c';'R11319_c';'R00264_c';'R03283_c';'R00026_c';'R03527_c';'R04949_c';'R04998_c';'R10035_c';'R10039_c';'R10040_c';'R02985_c';'R02558_c';'R02887_c';'R09375_c';'R09376_c';'R01177_c';'R01896_c';'R08240_c';'R04951_c';'R02569_c';'R00432_c';'R00727_c';'R03664_c';'R11896_c';'R02739_c';'R07346_c';'R02933_c';'R03751_c';'R00802_c';'R06088_c';'R11217_c';'R11219_c';'R04930_c';'R04770_c';'R09366_c';'R00782_c';'R02408_c';'R08193_c';'R00256_c';'R01395_c';'R00575_c';'R10948_c';'R10949_c';'R00717_c';'R01388_c';'R05922_c';'R08549_c';'R00548_c';'R05919_c';'R01000_c';'R03104_c';'R10563_c';'R10147_c';'R01512_c';'R07324_c';'R01993_c';'R03546_c';'R10079_c';'R09365_c';'R05086_c';'R10124_c';'R01078_c';'R01718_c';'R00801_c';'r726';'R06199_c';'R00985_c';'R00236_c';'R00316_c';'R00925_c';'R00926_c';'R01354_c';'R04880_c';'R05233_c';'R05234_c';'R06917_c';'R06927_c';'R07105_c';'R08281_c';'R08306_c';'R08310_c';'R02124_c';'R00754_c';'R00094_c';'R01776_c';'R01777_c';'R00425_c';'r503';'r504';'r505';'r506';'r1398';'r2083';'r564';'r353';'r1743';'r1903';'r1910';'r1917';'r1924';'r307';'r1518';'r1535';'r687';'r930';'r934';'r1941';'r680';'r676';'r1491';'r1667';'r996';'r855';'r193';'r351';'r718';'r570';'r573';'r576';'r579';'r868';'r527';'r1454';'r206';'r660';'r806';'r788';'r793';'r376';'r538';'r1664';'r1451';'r781';'r782';'r783';'r595';'r571';'r574';'r577';'r580';'r728';'r592';'r133';'r873';'r216';'r1523';'r812';'r817';'r2095';'r103';'r1780';'r623';'r112';'r771';'r91';'r878';'r889';'r1938';'r2161';'r90';'r989';'r140';'r141';'r501';'r394';'r1713';'r1831';'r1615';'r1541';'r1653';'r407';'r1860';'r533';'r1796';'R04360_c';'R00286_c';'R07676_c';'R06622_c';'R03566_c';'R06519_c';'R02204_c';'R11165_c';'R11220_c';'R01169_c';'R01989_c';'R04811_c';'R07769_c';'R09419_c';'R10587_c';'R06527_c';'R01461_c';'R06128_c';'R09676_c';'R06238_c';'R06261_c';'R10209_c';'R05976_c';'R07273_c';'R10951_c';'R01059_c';'R03819_c';'R10548_c';'R01547_c';'R04247_c';'R02590_c';'R01618_c';'R09726_c';'R03105_c';'R07461_c';'R07768_c';'R02687_c';'R03416_c';'R03417_c';'R11680_c';'R09827_c';'R03346_c';'R10952_c';'R02541_c';'R01470_c';'R06722_c';'R05981_c';'R06200_c';'R06516_c';'R01647_c';'R00470_c';'R10092_c';'R10827_c';'R01055_c';'R04773_c';'R01724_c';'R03905_c';'R08107_c';'R09073_c';'R10633_c';'r1568';'r978';'r981';'r1313';'r1328';'r2012';'r2064';'R10120_c';'R07763_c';'R10788_c';'R10995_c';'R10851_c';'R10852_c';'R09449_c';'R10828_c';'R03643_c';'R03816_c';'R04866_c';'R04867_c';'R07620_c';'R11399_c';'R01926_c';'R02052_c';'R07381_c';'R02051_c';'R07385_c';'R03980_c';'R04856_c';'R01421_c';'R01515_c';'R04804_c';'R07484_c';'R07771_c';'R11143_c';'R10822_c';'R10823_c';'r1844';'r1849';'r1854';'r1602';'r1616';'r1629';'r1983';'r1992';'r2018';'R07486_c';'R07491_c';'R07505_c';'R01456_c';'R07487_c';'R07492_c';'R02497_c';'R08954_c';'R10242_c';'R02661_c';'R00924_c';'R04432_c';'R02613_c';'R02382_c';'R04300_c';'R00277_c';'R00278_c';'R01710_c';'R04920_c';'R06364_c';'R07384_c';'R09413_c';'R09414_c';'R09415_c';'R06915_c';'R06936_c';'R06939_c';'R05632_c';'R00731_c';'R02363_c';'R04693_c';'R04884_c';'R05800_c';'R05801_c';'R10065_c';'R10953_c';'R11892_c';'R00512_c';'R01665_c';'R00158_c';'R04534_c';'R04536_c';'R04543_c';'R04566_c';'R04953_c';'R04258_c';'R05299_c';'R08114_c';'R08115_c';'R09134_c';'R09035_c';'R09036_c';'R01318_c';'R04480_c';'R03437_c';'R02352_c';'R02353_c';'R04681_c';'R04682_c';'R08945_c';'R08980_c';'R01278_c';'R03776_c';'R03856_c';'R03989_c';'R04753_c';'R06985_c';'R01175_c';'R01279_c';'R03990_c';'R03777_c';'R03857_c';'R04754_c';'R01317_c';'R02053_c';'R07064_c';'R07379_c';'R07387_c';'R07859_c';'r1645';'r1646';'r1647';'r1648';'r1649';'r1650';'R02920_c';'R03304_c';'R04301_c';'R04762_c';'R04764_c';'R04881_c';'R04887_c';'R05333_c';'R07493_c';'R11096_c';'R01457_c';'R03689_c';'R05703_c';'R07498_c';'R07499_c';'R07507_c';'r2003';'r2023';'r2027';'r2031';'r2040';'r2070';'r2052';'r2056';'r2060';'R05510_c';'R05511_c';'R06835_c';'R06838_c';'R08120_c';'R08121_c';'R09136_c';'R09220_c';'R09222_c';'R01103_c';'R01104_c';'R01194_c';'R01329_c';'R02926_c';'R03634_c';'R04019_c';'R04470_c';'R05549_c';'R06091_c';'R05961_c';'R02908_c';'R02919_c';'R04025_c';'R04890_c';'R04893_c';'R04894_c';'R04907_c';'R08346_c';'R08347_c';'R08348_c';'R02173_c';'R04674_c';'R04137_c';'R04204_c';'R04224_c';'R05595_c';'R06942_c';'R03045_c';'R02685_c';'R04170_c';'R04738_c';'R04740_c';'R04744_c';'R04746_c';'R04749_c';'R07002_c';'R07003_c';'R07004_c';'R07023_c';'R07024_c';'R07025_c';'R07026_c';'R07069_c';'R07070_c';'R07083_c';'R07084_c';'R07091_c';'R07092_c';'R07093_c';'R07094_c';'R07100_c';'R07113_c';'R07116_c';'R08280_c';'R09409_c';'R11905_c'};
reducedModel=removeReactions(refineModel,duplicateRxns,true,...
          true,true);
% unbound some reactions
irreversibleRxnsunbound = {'R00078_c';'R00100_c';'R00102_c';'R00126_c';'R00127_c';'R00132_c';'R00196_c';'R00282_c';'R00317_c';'R00381_c';'R00382_c';'R00471_c';'R00565_c';'R00597_c';'R00608_c';'R00610_c';'R00815_c';'R00817_c';'R00818_c';'R00821_c';'R00866_c';'R00966_c';'R01005_c';'R01025_c';'R01030_c';'R01056_c';'R01274_c';'R01280_c';'R01295_c';'R01312_c';'R01315_c';'R01352_c';'R01408_c';'R01416_c';'R01451_c';'R01462_c';'R01494_c';'R01514_c';'R01560_c';'R01626_c';'R01645_c';'R01708_c';'R01711_c';'R01827_c';'R01931_c';'R02025_c';'R02027_c';'R02057_c';'R02144_c';'R02208_c';'R02231_c';'R02239_c';'R02241_c';'R02250_c';'R02323_c';'R02383_c';'R02422_c';'R02480_c';'R02529_c';'R02532_c';'R02534_c';'R02556_c';'R02656_c';'R02746_c';'R02747_c';'R02756_c';'R02964_c';'R03026_c';'R03057_c';'R03172_c';'R03308_c';'R03315_c';'R03332_c';'R03353_c';'R03409_c';'R03435_c';'R03451_c';'R03478_c';'R03522_c';'R03659_c';'R03719_c';'R03772_c';'R03789_c';'R03805_c';'R03828_c';'R03875_c';'R03893_c';'R03938_c';'R03942_c';'R04072_c';'R04095_c';'R04212_c';'R04216_c';'R04404_c';'R04420_c';'R04452_c';'R04496_c';'R04533_c';'R04591_c';'R04620_c';'R04734_c';'R04751_c';'R04756_c';'R04858_c';'R04929_c';'R04964_c';'R05202_c';'R05287_c';'R05556_c';'R05616_c';'R05802_c';'R05916_c';'R05917_c';'R05918_c';'R05920_c';'R05923_c';'R05924_c';'R05970_c';'R05972_c';'R05973_c';'R05979_c';'R05980_c';'R05982_c';'R05986_c';'R06127_c';'R06171_c';'R06258_c';'R06259_c';'R06260_c';'R06262_c';'R06263_c';'R06264_c';'R06513_c';'R06517_c';'R06518_c';'R06601_c';'R07162_c';'R07215_c';'R07459_c';'R07488_c';'R07495_c';'R07758_c';'R07759_c';'R07760_c';'R07761_c';'R07762_c';'R07766_c';'R07767_c';'R07770_c';'R07816_c';'R08266_c';'R08555_c';'R08599_c';'R08701_c';'R09030_c';'R09034_c';'R09037_c';'R09038_c';'R09072_c';'R09087_c';'R09106_c';'R09289_c';'R09304_c';'R09395_c';'R09412_c';'R09597_c';'R09735_c';'R09782_c';'R09844_c';'R09845_c';'R09944_c';'R10116_c';'R10151_c';'R10231_c';'R10309_c';'R10455_c';'R10463_c';'R10532_c';'R10550_c';'R10565_c';'R10586_c';'R10648_c';'R10685_c';'R10686_c';'R10722_c';'R10815_c';'R10826_c';'R11014_c';'R11023_c';'R11166_c';'R11180_c';'R11216_c';'R11218_c';'R11221_c';'R11308_c';'R11583_c';'R11671_c';'R11861_c';'R11891_c';'R11929_c';'R12085_c';'r1113';'r1283';'r1487';'r1495';'r1510';'r1512';'r1526';'r1528';'r1536';'r1546';'r1559';'r1560';'r1561';'r1563';'r1564';'r1569';'r1570';'r1572';'r1576';'r1577';'r1584';'r1587';'r1588';'r1589';'r1644';'r1651';'r171';'r1718';'r1719';'r1727';'r1751';'r1756';'r1784';'r1786';'r1797';'r1800';'r1808';'r1809';'r1814';'r1827';'r1829';'r1839';'r1870';'r1877';'r1988';'r199';'r1997';'r2007';'r203';'r2035';'r204';'r2044';'r2074';'r2093';'r2099';'r2101';'r269';'r388';'r395';'r419';'r427';'r440';'r441';'r443';'r447';'r466';'r473';'r477';'r48';'r502';'r507';'r508';'r509';'r532';'r535';'r561';'r64';'r734';'r765';'r777';'r839';'r867';'r923';'r925';'r975';'r730';'R00533_c';'r1399';'r697';'r668';'R01009_c';'R00086_c';'r499'}; 
unboundModel=setParam(reducedModel,'ub',irreversibleRxnsunbound,1000);
unboundModel=setParam(unboundModel,'lb',irreversibleRxnsunbound,0);
model_iNR1320 = unboundModel;

%% STEP 4.2: MODEL REVISION 
% Updated by Nachon 29-OCT-2019
% prepare new reactions, genes, metabolites and EC numbers
load('ComplementaryData\iNR1320.mat');
rxnsToRemove = {'r0003';'r0004';'r0005';'r0006';'r0007';'r0008';'r0009';'r0010';'r0011';'r0012';'r0014';'r0015';'r0016';'r0017';'r0018';'r0019';'r0020';'r0022';'r0023';'r0024';'r0027';'r0028';'r0029';'r0030';'r0031';'r0032';'r0033';'r0034';'r0037';'r0038';'r0039';'r0041';'r0044';'r0046';'r0047';'r0048';'r0049';'r0050';'r0051';'r0052';'r0053';'r0054';'r0057';'r0058';'r0059';'r0060';'r0064';'r0065';'r0066';'r0067';'r0069';'r0070';'r0071';'r0072';'r0073';'r0074';'r0075';'r0076';'r0079';'r0080';'r0081';'r0082';'r0083';'r0085';'r0086';'r0087';'r0088';'r0089';'r0090';'r0091';'r0092';'r0093';'r0094';'r0095';'r0096';'r0097';'r0098';'r0099';'r0100';'r0101';'r0102';'r0103';'r0104';'r0105';'r0106';'r0108';'r0109';'r0110';'r0112';'r0113';'r0114';'r0115';'r0117';'r0119';'r0120';'r0122';'r0123';'r0124';'r0125';'r0126';'r0128';'r0129';'r0131';'r0135';'r0136';'r0137';'r0138';'r0139';'r0140';'r0141';'r0147';'r0149';'r0150';'r0151';'r0152';'r0155';'r0156';'r0157';'r0158';'r0159';'r0160';'r0162';'r0163';'r0164';'r0165';'r0166';'r0167';'r0168';'r0169';'r0170';'r0171';'r0172';'r0173';'r0176';'r0177';'r0178';'r0179';'r0180';'r0182';'r0183';'r0184';'r0187';'r0188';'r0189';'r0190';'r0192';'r0194';'r0195';'r0196';'r0197';'r0198';'r0199';'r0203';'r0204';'r0206';'r0207';'r0208';'r0209';'r0210';'r0211';'r0212';'r0215';'r0216';'r0217';'r0219';'r0224';'r0225';'r0226';'r0227';'r0228';'r0229';'r0230';'r0231';'r0233';'r0234';'r0235';'r0237';'r0238';'r0239';'r0240';'r0241';'r0242';'r0243';'r0245';'r0246';'r0247';'r0248';'r0251';'r0252';'r0255';'r0256';'r0257';'r0259';'r0260';'r0261';'r0262';'r0263';'r0264';'r0265';'r0266';'r0267';'r0268';'r0269';'r0270';'r0272';'r0273';'r0274';'r0275';'r0280';'r0281';'r0282';'r0283';'r0284';'r0287';'r0288';'r0289';'r0290';'r0291';'r0293';'r0294';'r0295';'r0296';'r0297';'r0298';'r0299';'r0301';'r0302';'r0303';'r0304';'r0305';'r0306';'r0307';'r0308';'r0309';'r0310';'r0311';'r0314';'r0315';'r0316';'r0319';'r0321';'r0322';'r0323';'r0324';'r0329';'r0330';'r0331';'r0334';'r0336';'r0337';'r0338';'r0339';'r0340';'r0341';'r0342';'r0343';'r0344';'r0345';'r0346';'r0347';'r0348';'r0349';'r0350';'r0351';'r0352';'r0353';'r0354';'r0355';'r0356';'r0357';'r0358';'r0359';'r0360';'r0361';'r0362';'r0372';'r0373';'r0374';'r0375';'r0376';'r0377';'r0378';'r0379';'r0381';'r0382';'r0383';'r0384';'r0385';'r0386';'r0387';'r0388';'r0389';'r0390';'r0391';'r0392';'r0393';'r0394';'r0395';'r0397';'r0398';'r0399';'r0400';'r0401';'r0402';'r0403';'r0408';'r0409';'r0410';'r0411';'r0412';'r0414';'r0415';'r0416';'r0417';'r0418';'r0420';'r0421';'r0422';'r0423';'r0424';'r0425';'r0426';'r0427';'r0428';'r0429';'r0430';'r0431';'r0432';'r0433';'r0434';'r0435';'r0436';'r0437';'r0438';'r0439';'r0441';'r0442';'r0443';'r0444';'r0445';'r0446';'r0447';'r0448';'r0450';'r0451';'r0452';'r0453';'r0454';'r0455';'r0456';'r0457';'r0458';'r0461';'r0463';'r0464';'r0465';'r0470';'r0471';'r0472';'r0473';'r0474';'r0475';'r0476';'r0478';'r0479';'r0480';'r0481';'r0482';'r0483';'r0485';'r0486';'r0487';'r0488';'r0491';'r0496';'r0497';'r0498';'r0499';'r0500';'r0501';'r0502';'r0503';'r0506';'r0510';'r0511';'r0512';'r0514';'r0515';'r0516';'r0517';'r0518';'r0519';'r0520';'r0521';'r0525';'r0526';'r0527';'r0528';'r0529';'r0531';'r0532';'r0533';'r0535';'r0536';'r0537';'r0538';'r0539';'r0540';'r0541';'r0542';'r0543';'r0544';'r0545';'r0546';'r0547';'r0548';'r0549';'r0550';'r0551';'r0552';'r0554';'r0555';'r0557';'r0558';'r0559';'r0560';'r0561';'r0562';'r0563';'r0564';'r0565';'r0566';'r0567';'r0568';'r0569';'r0570';'r0571';'r0572';'r0573';'r0574';'r0575';'r0576';'r0577';'r0578';'r0579';'r0580';'r0581';'r0583';'r0584';'r0585';'r0586';'r0590';'r0591';'r0592';'r0593';'r0594';'r0595';'r0596';'r0597';'r0598';'r0599';'r0600';'r0601';'r0602';'r0603';'r0604';'r0605';'r0606';'r0607';'r0608';'r0610';'r0611';'r0612';'r0613';'r0614';'r0617';'r0618';'r0619';'r0620';'r0621';'r0622';'r0623';'r0624';'r0625';'r0627';'r0628';'r0629';'r0630';'r0631';'r0632';'r0633';'r0634';'r0635';'r0636';'r0638';'r0639';'r0640';'r0641';'r0642';'r0643';'r0644';'r0645';'r0647';'r0648';'r0650';'r0651';'r0652';'r0653';'r0654';'r0655';'r0656';'r0657';'r0658';'r0659';'r0660';'r0661';'r0662';'r0663';'r0664';'r0665';'r0666';'r0667';'r0669';'r0670';'r0671';'r0672';'r0673';'r0674';'r0675';'r0676';'r0677';'r0678';'r0679';'r0680';'r0681';'r0682';'r0684';'r0685';'r0688';'r0689';'r0690';'r0691';'r0693';'r0694';'r0695';'r0696';'r0697';'r0698';'r0699';'r0700';'r0702';'r0705';'r0706';'r0707';'r0708';'r0709';'r0710';'r0711';'r0712';'r0713';'r0714';'r0715';'r0717';'r0718';'r0719';'r0720';'r0721';'r0722';'r0724';'r0725';'r0726';'r0727';'r0728';'r0729';'r0730';'r0731';'r0732';'r0733';'r0734';'r0735';'r0736';'r0737';'r0738';'r0739';'r0740';'r0743';'r0744';'r0746';'r0747';'r0748';'r0749';'r0750';'r0752';'r0753';'r0755';'r0756';'r0757';'r0758';'r0759';'r0762';'r0763';'r0764';'r0765';'r0766';'r0767';'r0768';'r0769';'r0770';'r0771';'r0772';'r0773';'r0774';'r0775';'r0776';'r0777';'r0778';'r0779';'r0780';'r0781';'r0782';'r0783';'r0784';'r0785';'r0786';'r0787';'r0788';'r0789';'r0790';'r0791';'r0792';'r0793';'r0794';'r0795';'r0796';'r0797';'r0798';'r0799';'r0808';'r0821';'r0822';'r0823';'r0824';'r0825';'r0826';'r0827';'r0828';'r0829';'r0830';'r0831';'r0832';'r0833';'r0834';'r0835';'r0836';'r0837';'r0838';'r0839';'r0840';'r0841';'r0842';'r0843';'r0844';'r0845';'r0846';'r0847';'r0848';'r0849';'r0850';'r0851';'r0852';'r0853';'r0854';'r0855';'r0856';'r0857';'r0858';'r0859';'r0860';'r0861';'r0862';'r0863';'r0864';'r0865';'r0866';'r0867';'r0868';'r0869';'r0870';'r0871';'r0872';'r0873';'r0874';'r0875';'r0876';'r0877';'r0878';'r0879';'r0880';'r0881';'r0882';'r0883';'r0884';'r0885';'r0886';'r0887';'r0888';'r0889';'r0890';'r0891';'r0892';'r0893';'r0894';'r0895';'r0896';'r0897';'r0898';'r0899';'r0900';'r0901';'r0902';'r0903';'r0904';'r0905';'r0906';'r0907';'r0908';'r0909';'r0910';'r0911';'r0912';'r0913';'r0914';'r0915';'r0916';'r0917';'r0918';'r0919';'r0920';'r0921';'r0922';'r0923';'r0924';'r0925';'r0926';'r0927';'r0928';'r0929';'r0930';'r0931';'r0932';'r0933';'r0934';'r0935';'r0936';'r0937';'r0938';'r0939';'r0940';'r0941';'r0942';'r0943';'r0944';'r0945';'r0946';'r0947';'r0948';'r0949';'r0950';'r0951';'r0952';'r0953';'r0954';'r0955';'r0956';'r0957';'r0958';'r0960';'r0961';'r0962';'r0964';'r0965';'r0966';'r0968';'r0969';'r0970';'r0972';'r0973';'r0974';'r0975';'r0976';'r0979';'r0980';'r0981';'r0982';'r0983';'r0984';'r0986';'r0987';'r0988';'r0990';'r0991';'r0992';'r0994';'r0995';'r0996';'r0998';'r0999';'r1000';'r1002';'r1003';'r1004';'r1005';'r1006';'r1007';'r1008';'r1009';'r1010';'r1011';'r1012';'r1013';'r1014';'r1015';'r1016';'r1017';'r1018';'r1019';'r1020';'r1021';'r1022';'r1023';'r1024';'r1025';'r1026';'r1027';'r1028';'r1029';'r1030';'r1031';'r1032';'r1033';'r1034';'r1035';'r1036';'r1037';'r1038';'r1039';'r1040';'r1041';'r1042';'r1043';'r1044';'r1045';'r1046';'r1047';'r1048';'r1049';'r1050';'r1051';'r1052';'r1053';'r1054';'r1055';'r1056';'r1057';'r1058';'r1059';'r1060';'r1061';'r1062';'r1063';'r1064';'r1065';'r1066';'r1067';'r1068';'r1069';'r1070';'r1071';'r1072';'r1073';'r1074';'r1075';'r1076';'r1077';'r1078';'r1079';'r1080';'r1081';'r1082';'r1083';'r1084';'r1085';'r1086';'r1087';'r1088';'r1089';'r1090';'r1091';'r1092';'r1093';'r1094';'r1095';'r1096';'r1097';'r1098';'r1099';'r1100';'r1101';'r1102';'r1103';'r1104';'r1105';'r1106';'r1107';'r1108';'r1109';'r1110';'r1111';'r1112';'r1114';'r1115';'r1116';'r1117';'r1118';'r1119';'r1120';'r1121';'r1122';'r1123';'r1125';'r1126';'r1127';'r1128';'r1129';'r1131';'r1132';'r1133';'r1134';'r1135';'r1136';'r1137';'r1138';'r1139';'r1140';'r1141';'r1142';'r1143';'r1144';'r1145';'r1146';'r1147';'r1149';'r1150';'r1152';'r1154';'r1158';'r1160';'r1161';'r1162';'r1177';'r1178';'r1179';'r1182';'r1183';'r1185';'r1188';'r1195';'r1196';'r1197';'r1198';'r1205';'r1217';'r1218';'r1220';'r1223';'r1226';'r1227';'r1228';'r1229';'r1230';'r1231';'r1232';'r1234';'r1237';'r1239';'r1241';'r1242';'r1262';'r1265';'r1266';'r1267';'r1269';'r1274';'r1276';'r1277';'r1278';'r1280';'r1282';'r1285';'r1287';'r1288';'r1315';'r1316';'r1317';'r1322';'r1324';'r1332';'r1333';'r1336';'r1345';'r1356';'r1359';'r1360';'r1361';'r1363';'r1366';'r1367';'r1368';'r1371';'r1372';'r1373';'r1374';'r1376';'r1377';'r1379';'r1381';'r1382';'r1383';'r1384';'r1386';'r1391';'r1394';'r1410';'r1411';'r1413';'r0703';'r0704';'r0084';'r0687';'r0146';'r0686';'r0127';'r1148';'r0193';'r0025';'r0553';'r0062';'r0063';'r0467';'r0462';'r0143';'r0144';'r0468';'r0469';'r0111';'r0368';'r1124';'r0523';'r0413';'r0419';'r0683';'r1155';'r0134';'r0477';'r0145';'r0701';'r0716';'r0013';'r0186';'r0153';'r0142';'r0524';'r0118';'r0244';'r0508';'r0513';'r0258';'r0300';'r0466';'r0609';'r0646';'r0522';'r0490';'r0587';'r0250';'r0277';'r0249';'biomass';'protein';'dna';'rna';'lipid';'carbohydrate';'matp';'3amp';'r0001';'13glucanIN';'4hpoaIN';'acIN';'adIN';'akgIN';'alaIN';'amaIN';'anIN';'aolIN';'arabIN';'arabinIN';'argIN';'asnIN';'aspIN';'balaIN';'bdglcIN';'bmOUT';'c100IN';'c120IN';'c140IN';'c150IN';'c160IN';'c161IN';'c162IN';'c170IN';'c171IN';'c180IN';'c181IN';'c182IN';'c183IN';'c200IN';'c40IN';'c50IN';'c60IN';'c70IN';'c80IN';'c90IN';'cb15lctIN';'cellobIN';'celluIN';'chibIN';'chitIN';'chitoIN';'choIN';'citrIN';'co2OUT';'cordycepinOUT';'cyneIN';'cysIN';'cytsIN';'dglcIN';'ethnitIN';'fmnIN';'forIN';'fruIN';'fumIN';'gabaIN';'glacIN';'glcIN';'glcn15lacIN';'glcnIN';'glIN';'glnIN';'gluIN';'glyaIN';'glycogenIN';'glyIN';'gnIN';'h2oOUT';'h2sIN';'h2so3IN';'hcysIN';'hisIN';'hno2IN';'hno3IN';'hyxnIN';'icitIN';'idolIN';'ileIN';'lactIN';'laolIN';'larabIN';'leuIN';'llacIN';'lrlIN';'lysIN';'malIN';'manIN';'mannanIN';'meliIN';'metholIN';'metIN';'mltIN';'mltioseIN';'mntIN';'myoiIN';'nagIN';'nh3IN';'nicaIN';'nicdIN';'nr0010';'o2IN';'oaIN';'ornIN';'oxalIN';'paaIN';'pabaIN';'pheIN';'piIN';'pimIN';'poaIN';'proIN';'propIN';'pyrIN';'quinIN';'r0002';'r0021';'r0130';'r0133';'r0175';'r0181';'r0191';'r0205';'r0223';'r0312';'r0460';'r0489';'r0492';'r0505';'r0507';'r0509';'r0589';'r0649';'r1164';'r1165';'r1166';'r1167';'r1168';'r1169';'r1170';'r1171';'r1172';'r1173';'r1174';'r1175';'r1176';'r1180';'r1181';'r1184';'r1186';'r1187';'r1189';'r1190';'r1191';'r1192';'r1193';'r1194';'r1199';'r1200';'r1201';'r1202';'r1203';'r1204';'r1206';'r1207';'r1208';'r1209';'r1210';'r1211';'r1212';'r1213';'r1214';'r1215';'r1216';'r1219';'r1221';'r1222';'r1224';'r1225';'r1233';'r1235';'r1236';'r1238';'r1240';'r1243';'r1244';'r1245';'r1246';'r1247';'r1248';'r1249';'r1250';'r1251';'r1252';'r1253';'r1258';'r1259';'r1260';'r1263';'r1264';'r1270';'r1272';'r1273';'r1275';'r1284';'r1286';'r1289';'r1290';'r1291';'r1292';'r1293';'r1294';'r1295';'r1296';'r1297';'r1298';'r1299';'r1300';'r1301';'r1302';'r1303';'r1304';'r1305';'r1308';'r1309';'r1314';'r1318';'r1321';'r1323';'r1325';'r1326';'r1327';'r1329';'r1330';'r1331';'r1334';'r1335';'r1337';'r1338';'r1339';'r1340';'r1341';'r1342';'r1343';'r1344';'r1346';'r1347';'r1348';'r1349';'r1350';'r1352';'r1353';'r1355';'r1357';'r1358';'r1362';'r1364';'r1365';'r1369';'r1370';'r1375';'r1378';'r1380';'r1385';'r1387';'r1388';'r1389';'r1390';'r1392';'r1393';'r1396';'r1397';'r1400';'r1402';'r1403';'r1404';'r1405';'r1406';'r1407';'r1408';'r1409';'r1412';'r1414';'r1415';'r1416';'r1417';'r1418';'r1419';'r1420';'r1421';'r1422';'r1423';'r1424';'r1425';'r1427';'r1428';'r1429';'r1430';'r1431';'r1432';'r1433';'r1434';'r1435';'r1436';'r1437';'r1438';'r1439';'r1440';'r1441';'r1442';'r1443';'r1444';'r1445';'r1446';'r1447';'r1448';'r1449';'r1450';'raffIN';'rgtIN';'ribIN';'serIN';'slfIN';'sorIN';'starIN';'succIN';'sucIN';'sulfurIN';'thmIN';'thrIN';'treIN';'trpIN';'tyrIN';'uraIN';'ureaIN';'uriIN';'valIN';'xanIN';'xolIN';'xylanIN';'xylIN';'r48';'r64';'r171';'r199';'r203';'r204';'r269';'r388';'r395';'r419';'r427';'r440';'r441';'r443';'r447';'r466';'r473';'r477';'r499';'r502';'r507';'r508';'r509';'r532';'r535';'r561';'r668';'r697';'r730';'r734';'r765';'r777';'r839';'r867';'r923';'r925';'r975';'r1113';'r1283';'r1399';'r1487';'r1495';'r1510';'r1512';'r1526';'r1528';'r1536';'r1540';'r1542';'r1546';'r1559';'r1560';'r1561';'r1563';'r1564';'r1569';'r1570';'r1572';'r1576';'r1577';'r1584';'r1587';'r1588';'r1589';'r1644';'r1651';'r1718';'r1719';'r1727';'r1751';'r1756';'r1784';'r1786';'r1797';'r1800';'r1808';'r1809';'r1814';'r1827';'r1829';'r1839';'r1870';'r1877';'r1988';'r1997';'r2007';'r2035';'r2044';'r2074';'r2093';'r2099';'r2101';'R00078_c';'R00086_c';'R00100_c';'R00102_c';'R00126_c';'R00127_c';'R00132_c';'R00196_c';'R00282_c';'R00317_c';'R00351_c';'R00471_c';'R00533_c';'R00565_c';'R00597_c';'R00608_c';'R00610_c';'R00815_c';'R00817_c';'R00818_c';'R00821_c';'R00866_c';'R00966_c';'R01005_c';'R01009_c';'R01025_c';'R01030_c';'R01056_c';'R01101_c';'R01274_c';'R01280_c';'R01295_c';'R01312_c';'R01315_c';'R01352_c';'R01408_c';'R01416_c';'R01451_c';'R01462_c';'R01494_c';'R01514_c';'R01560_c';'R01626_c';'R01645_c';'R01663_c';'R01708_c';'R01711_c';'R01827_c';'R01845_c';'R01931_c';'R02025_c';'R02027_c';'R02057_c';'R02144_c';'R02208_c';'R02231_c';'R02239_c';'R02241_c';'R02250_c';'R02323_c';'R02383_c';'R02422_c';'R02480_c';'R02529_c';'R02532_c';'R02534_c';'R02556_c';'R02656_c';'R02746_c';'R02747_c';'R02756_c';'R02964_c';'R03026_c';'R03057_c';'R03172_c';'R03308_c';'R03315_c';'R03332_c';'R03353_c';'R03409_c';'R03435_c';'R03451_c';'R03478_c';'R03522_c';'R03659_c';'R03719_c';'R03772_c';'R03789_c';'R03805_c';'R03828_c';'R03875_c';'R03893_c';'R03938_c';'R03942_c';'R04072_c';'R04095_c';'R04212_c';'R04216_c';'R04404_c';'R04420_c';'R04452_c';'R04496_c';'R04533_c';'R04591_c';'R04620_c';'R04734_c';'R04751_c';'R04756_c';'R04858_c';'R04929_c';'R04964_c';'R05202_c';'R05287_c';'R05616_c';'R05802_c';'R05916_c';'R05917_c';'R05918_c';'R05920_c';'R05923_c';'R05924_c';'R05970_c';'R05972_c';'R05973_c';'R05979_c';'R05980_c';'R05982_c';'R05986_c';'R06127_c';'R06171_c';'R06258_c';'R06259_c';'R06260_c';'R06262_c';'R06263_c';'R06264_c';'R06513_c';'R06517_c';'R06518_c';'R06601_c';'R07162_c';'R07215_c';'R07459_c';'R07488_c';'R07495_c';'R07758_c';'R07759_c';'R07760_c';'R07761_c';'R07762_c';'R07766_c';'R07767_c';'R07770_c';'R07816_c';'R08266_c';'R08555_c';'R08599_c';'R09030_c';'R09034_c';'R09037_c';'R09038_c';'R09072_c';'R09087_c';'R09106_c';'R09289_c';'R09304_c';'R09395_c';'R09412_c';'R09597_c';'R09735_c';'R09782_c';'R09844_c';'R09845_c';'R09944_c';'R10116_c';'R10151_c';'R10231_c';'R10309_c';'R10455_c';'R10463_c';'R10532_c';'R10550_c';'R10565_c';'R10586_c';'R10648_c';'R10685_c';'R10722_c';'R10815_c';'R10826_c';'R11014_c';'R11023_c';'R11166_c';'R11180_c';'R11216_c';'R11218_c';'R11221_c';'R11308_c';'R11583_c';'R11671_c';'R11861_c';'R11891_c';'R11929_c';'R12085_c';'R09144_c';'R00381_c';'R00382_c';'R05556_c';'R08701_c';'R10686_c'};
reducedModel = removeReactions(model,rxnsToRemove,true,...
            false,false);

[~, newMets]=xlsread('ComplementaryData\plus1329.xlsx','mets');
metsToAdd = struct();
metsToAdd.mets = newMets(2:end,1);
metsToAdd.metNames = newMets(2:end,2);
metsToAdd.compartments = 'c';
metPlusModel=addMets(reducedModel,metsToAdd);

[~, SheetS]=xlsread('ComplementaryData\plus1329.xlsx','Sheet1');
rxnsToAdd = struct();
rxnsToAdd.rxns = SheetS(2:end,1);
rxnsToAdd.rxnNames = SheetS(2:end,2);
rxnsToAdd.equations = SheetS(2:end,3);
rxnsToAdd.grRules = SheetS(2:end,5);
rxnsToAdd.eccodes = SheetS(2:end,4);
rxnsToAdd.rxnNotes = SheetS(2:end,6);

newModel=addRxns(metPlusModel,rxnsToAdd,1,'c',true,true);
newModel.mets = newModel.metNames;
newModel.equations = constructEquations(newModel);
newRxnWv = newModel.rxns;

% prepare template model from model_iNR1320
badrxnsToRemove = {'R08701_c','R01009_c','R11308_c','r1487','r1495','R02756_c','R03875_c','R05916_c','R05917_c','R05918_c','R05920_c','R05923_c','R05924_c','R05970_c','R05972_c','R05973_c','R05979_c','R05982_c','R05986_c','R06258_c','R06259_c','R06260_c','R07488_c','R07758_c','R07759_c','R07760_c','R07762_c','R07770_c','R07816_c','R08599_c','R09037_c','R09304_c','R09412_c','R10151_c','R11671_c','R10686_c'};
goodModel = removeReactions(model_iNR1320,badrxnsToRemove,true,...
            false,true);

% add reactions, genes and metabolites into the prepared template
kmodel=addRxnsGenesMets(goodModel,newModel,newRxnWv,true,true,4);
kmodel = setParam(kmodel,'ub',newRxnWv,1000);
kmodel = setParam(kmodel,'lb',newRxnWv,0);
%kmodel = setParam(kmodel,'ub','wv_R00057',0);

pmodel = contractModel(kmodel);
[a, b]=ismember(pmodel.rxns,rxnsToAdd.rxns);
I=find(a);

pmodel.rxnNotes(I)=rxnsToAdd.rxnNotes(b(I));
ppModel = contractModel(pmodel);

% recheck and remove unused genes
finalModel=deleteUnusedGenes(ppModel);

% recheck and remove some duplicated metabolites
[~, textData1]=xlsread('ComplementaryData\supplementary.xlsx','Table S2');
metNames = struct();
metNames.old = textData1(4:end,2);
metNames.new = textData1(4:end,3);
[a, b]=ismember(finalModel.metNames,metNames.old);
I=find(a);
finalModel.metNames(I)=metNames.new(b(I));
finalModel.metNames = lower(finalModel.metNames);
lastModel=contractModel(finalModel);


%% STEP 5: MODEL VALIDATION
% The model validation was performed against experimentation.
finalValidateModel = lastModel;
finalValidateModel = setParam(finalValidateModel,'lb',{'bmOUT','cordycepinOUT'},0);
finalValidateModel = setParam(finalValidateModel,'eq',{'matp'},1);
finalValidateModel = setParam(finalValidateModel,'ub',{'bmOUT','cordycepinOUT'},1000);
finalValidateModel = setParam(finalValidateModel,'obj',{'bmOUT'},1);
sol = solveLP(finalValidateModel,1);
% The following information was reported in Table 2
sugars = {'glcIN' 'fruIN' 'arabIN' 'xylIN' 'sucIN'};
uptake = [0.1448,0.1472,0.1074,0.0681,0.0815]; %observed from experiment
sumax(1) = sol;

for i = 1:numel(uptake)
    model=setParam(finalValidateModel,'ub',sugars,0);
    model=setParam(model,'ub',sugars(i),uptake(i));
    sugar = sugars(i);
    sol = solveLP(model,1);
    sumax(i) = sol;
    umax = num2str(sol.f*-1);
    fprintf([num2str(sol.f*-1) '\n']);
end

iNR1329 = contractModel(finalValidateModel)




% If the model prediction matched the experimentation with the Error rates below
% 5%, then we will export the model for the latter analysis.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iNR1329.id = 'iNR1329';
iNR1329.description = 'C. militaris GEM by Nachon Raethong [29-OCT-2019]';
iNR1329.annotation = [];
exportForGit(iNR1329,'model','C:\Users\Dell\Documents\GitHub\Cordyceps_militaris-GEM\');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

