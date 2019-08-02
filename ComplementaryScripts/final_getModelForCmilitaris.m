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

%% STEP 4: MODEL REFINEMENT
% Remove duplicated reactions
refineModel = modelBigger;
duplicateRxns = {'r12';'r13';'r120';'r135';'r158';'r185';'r186';'r188';'r191';'r198';'r236';'r238';'r250';'r274';'r334';'r372';'r374';'r377';'r378';'r379';'r458';'r464';'r495';'r517';'r531';'r559';'r609';'r647';'r666';'r685';'r686';'r714';'r717';'r729';'r749';'r754';'r769';'r786';'r787';'r798';'r799';'r800';'r801';'r803';'r804';'r805';'r807';'r811';'r813';'r816';'r819';'r828';'r830';'r831';'r833';'r835';'r836';'r837';'r845';'r857';'r871';'r888';'r896';'r897';'r899';'r900';'r924';'r929';'r939';'r947';'r958';'r964';'r966';'r969';'r971';'r977';'r1156';'r1157';'r1159';'r1163';'r1279';'r1281';'r1395';'r1459';'r1461';'r1464';'r1466';'r1496';'r1498';'r1500';'r1504';'r1507';'r1509';'r1515';'r1533';'r1575';'r1578';'r1581';'r1582';'r1583';'r1676';'r1685';'r1686';'r1688';'r1689';'r1693';'r1728';'r1734';'r1737';'r1739';'r1749';'r1762';'r1768';'r1769';'r1770';'r1779';'r1792';'r1798';'r1817';'r2098';'r2163';'r2179';'R00009_c';'R00028_c';'R00093_c';'R00103_c';'R00112_c';'R00114_c';'R00139_c';'R00156_c';'R00190_c';'R00197_c';'R00200_c';'R00238_c';'R00239_c';'R00245_c';'R00258_c';'R00268_c';'R00289_c';'R00291_c';'R00299_c';'R00310_c';'R00319_c';'R00330_c';'R00342_c';'R00355_c';'R00369_c';'R00372_c';'R00416_c';'R00426_c';'R00428_c';'R00430_c';'R00472_c';'R00479_c';'R00489_c';'R00513_c';'R00516_c';'R00517_c';'R00570_c';'R00588_c';'R00621_c';'R00658_c';'R00678_c';'R00691_c';'R00699_c';'R00702_c';'R00705_c';'R00706_c';'R00707_c';'R00708_c';'R00714_c';'R00715_c';'R00720_c';'R00722_c';'R00736_c';'R00742_c';'R00746_c';'R00755_c';'R00756_c';'R00760_c';'R00762_c';'R00765_c';'R00768_c';'R00774_c';'R00835_c';'R00842_c';'R00848_c';'R00851_c';'R00885_c';'R00891_c';'R00893_c';'R00896_c';'R00908_c';'R00922_c';'R00927_c';'R00935_c';'R00936_c';'R00937_c';'R00939_c';'R00940_c';'R00942_c';'R00945_c';'R00955_c';'R00959_c';'R00962_c';'R00967_c';'R00968_c';'R00970_c';'R00994_c';'R01004_c';'R01015_c';'R01041_c';'R01049_c';'R01057_c';'R01061_c';'R01067_c';'R01071_c';'R01072_c';'R01073_c';'R01082_c';'R01083_c';'R01086_c';'R01098_c';'R01123_c';'R01127_c';'R01135_c';'R01137_c';'R01138_c';'R01163_c';'R01184_c';'R01213_c';'R01224_c';'R01226_c';'R01229_c';'R01230_c';'R01248_c';'R01251_c';'R01252_c';'R01273_c';'R01324_c';'R01325_c';'R01351_c';'R01360_c';'R01364_c';'R01372_c';'R01398_c';'R01404_c';'R01481_c';'R01497_c';'R01518_c';'R01529_c';'R01542_c';'R01548_c';'R01549_c';'R01602_c';'R01641_c';'R01648_c';'R01678_c';'R01698_c';'R01701_c';'R01702_c';'R01752_c';'R01768_c';'R01769_c';'R01773_c';'R01775_c';'R01819_c';'R01830_c';'R01855_c';'R01857_c';'R01858_c';'R01868_c';'R01872_c';'R01880_c';'R01899_c';'R01900_c';'R01914_c';'R01933_c';'R01934_c';'R01936_c';'R01939_c';'R01940_c';'R01954_c';'R01960_c';'R01976_c';'R01978_c';'R01986_c';'R02016_c';'R02017_c';'R02019_c';'R02021_c';'R02024_c';'R02080_c';'R02082_c';'R02085_c';'R02091_c';'R02093_c';'R02094_c';'R02096_c';'R02097_c';'R02098_c';'R02103_c';'R02106_c';'R02107_c';'R02109_c';'R02112_c';'R02133_c';'R02222_c';'R02235_c';'R02236_c';'R02251_c';'R02282_c';'R02283_c';'R02291_c';'R02300_c';'R02301_c';'R02315_c';'R02320_c';'R02326_c';'R02327_c';'R02331_c';'R02332_c';'R02366_c';'R02367_c';'R02371_c';'R02372_c';'R02401_c';'R02410_c';'R02413_c';'R02425_c';'R02433_c';'R02466_c';'R02485_c';'R02487_c';'R02488_c';'R02508_c';'R02519_c';'R02521_c';'R02530_c';'R02536_c';'R02537_c';'R02566_c';'R02570_c';'R02571_c';'R02596_c';'R02619_c';'R02649_c';'R02662_c';'R02665_c';'R02668_c';'R02678_c';'R02695_c';'R02697_c';'R02701_c';'R02703_c';'R02720_c';'R02736_c';'R02740_c';'R02750_c';'R02864_c';'R02872_c';'R02874_c';'R02940_c';'R02957_c';'R02971_c';'R02978_c';'R02984_c';'R03004_c';'R03012_c';'R03024_c';'R03038_c';'R03051_c';'R03139_c';'R03158_c';'R03174_c';'R03180_c';'R03220_c';'R03222_c';'R03236_c';'R03237_c';'R03238_c';'R03239_c';'R03291_c';'R03293_c';'R03300_c';'R03302_c';'R03313_c';'R03316_c';'R03321_c';'R03348_c';'R03355_c';'R03443_c';'R03471_c';'R03530_c';'R03531_c';'R03596_c';'R03635_c';'R03646_c';'R03650_c';'R03652_c';'R03654_c';'R03655_c';'R03657_c';'R03658_c';'R03662_c';'R03663_c';'R03665_c';'R03778_c';'R03815_c';'R03858_c';'R03869_c';'R03909_c';'R03919_c';'R03920_c';'R03921_c';'R03936_c';'R03947_c';'R03968_c';'R03991_c';'R04001_c';'R04007_c';'R04027_c';'R04065_c';'R04097_c';'R04125_c';'R04138_c';'R04144_c';'R04171_c';'R04173_c';'R04178_c';'R04188_c';'R04209_c';'R04225_c';'R04355_c';'R04371_c';'R04378_c';'R04385_c';'R04386_c';'R04391_c';'R04424_c';'R04426_c';'R04428_c';'R04430_c';'R04439_c';'R04440_c';'R04444_c';'R04445_c';'R04506_c';'R04535_c';'R04537_c';'R04544_c';'R04559_c';'R04560_c';'R04568_c';'R04725_c';'R04726_c';'R04742_c';'R04747_c';'R04783_c';'R04862_c';'R04882_c';'R04883_c';'R04888_c';'R04889_c';'R04891_c';'R04892_c';'R04903_c';'R04909_c';'R04928_c';'R04936_c';'R04942_c';'R04944_c';'R04945_c';'R04946_c';'R04952_c';'R04954_c';'R04956_c';'R04957_c';'R04959_c';'R04960_c';'R04962_c';'R04963_c';'R04965_c';'R04967_c';'R04968_c';'R04970_c';'R04972_c';'R04982_c';'R04983_c';'R04988_c';'R04996_c';'R05046_c';'R05050_c';'R05051_c';'R05052_c';'R05064_c';'R05066_c';'R05068_c';'R05069_c';'R05071_c';'R05085_c';'R05231_c';'R05237_c';'R05238_c';'R05286_c';'R05487_c';'R05551_c';'R05576_c';'R05590_c';'R05614_c';'R05639_c';'R05640_c';'R05731_c';'R06004_c';'R06087_c';'R06114_c';'R06154_c';'R06366_c';'R06525_c';'R06526_c';'R06590_c';'R06740_c';'R06941_c';'R07034_c';'R07035_c';'R07104_c';'R07136_c';'R07168_c';'R07363_c';'R07364_c';'R07392_c';'R07411_c';'R07412_c';'R07443_c';'R07460_c';'R07481_c';'R07483_c';'R07494_c';'R07497_c';'R07599_c';'R07600_c';'R07601_c';'R07602_c';'R07603_c';'R07604_c';'R07618_c';'R07888_c';'R07891_c';'R07892_c';'R07895_c';'R07896_c';'R07899_c';'R07934_c';'R07937_c';'R07942_c';'R07950_c';'R07953_c';'R07977_c';'R07978_c';'R07979_c';'R07981_c';'R08090_c';'R08146_c';'R08218_c';'R08221_c';'R08231_c';'R08232_c';'R08235_c';'R08243_c';'R08244_c';'R08282_c';'R08283_c';'R08307_c';'R08363_c';'R08364_c';'R08381_c';'R08639_c';'R08734_c';'R08739_c';'R09076_c';'R09099_c';'R09245_c';'R09246_c';'R09300_c';'R09301_c';'R09372_c';'R09394_c';'R09783_c';'R09951_c';'R09993_c';'R10046_c';'R10052_c';'R10074_c';'R10115_c';'R10119_c';'R10170_c';'R10221_c';'R10619_c';'R10700_c';'R10712_c';'R10907_c';'R10996_c';'R10997_c';'R10998_c';'R11104_c';'R11262_c';'R11316_c';'R11329_c';'R11528_c';'R11529_c';'R11765_c';'R11783_c';'R11893_c';'R11894_c';'R11895_c';'R11906_c';'R02317_c';'R04639_c';'R05048_c';'R06982_c';'R08637_c';'R08698_c';'R11098_c';'R11099_c';'R04241_c';'R00698_c';'R04546_c';'R02918_c';'R00996_c';'R00332_c';'R02090_c';'R05921_c';'R04984_c';'R05641_c';'R10993_c';'R00348_c';'R00269_c';'R08359_c';'R01092_c';'R00006_c';'R00226_c';'R03050_c';'R04672_c';'R04673_c';'R08648_c';'R01068_c';'R01070_c';'R01829_c';'R02568_c';'R00709_c';'R09107_c';'R03661_c';'R03194_c';'R03648_c';'R01699_c';'R03270_c';'R00209_c';'R00014_c';'R05577_c';'R00137_c';'R02058_c';'R01800_c';'R00410_c';'R00248_c';'R00590_c';'R03601_c';'R04859_c';'R11593_c';'R10994_c';'R00830_c';'R00527_c';'R01818_c';'R03444_c';'R02472_c';'R01039_c';'R01799_c';'R05578_c';'R01262_c';'R01687_c';'R03867_c';'R03916_c';'R03970_c';'R03971_c';'R04935_c';'R00508_c';'R00177_c';'R04771_c';'R02164_c';'R02670_c';'R00602_c';'R02517_c';'R07129_c';'R02702_c';'R02909_c';'R03628_c';'R02396_c';'R00259_c';'R04941_c';'R00794_c';'R00796_c';'R07506_c';'R01186_c';'R03629_c';'R04121_c';'R05259_c';'R03098_c';'R04390_c';'R04863_c';'R01209_c';'R04441_c';'R05070_c';'R02111_c';'R03243_c';'R00694_c';'R00734_c';'R00187_c';'R01870_c';'R03660_c';'R00858_c';'R05963_c';'R07809_c';'R07810_c';'R10831_c';'R02199_c';'R10991_c';'R01214_c';'R01090_c';'R03656_c';'R03425_c';'R03033_c';'R11319_c';'R00264_c';'R03283_c';'R00026_c';'R03527_c';'R04949_c';'R04998_c';'R10035_c';'R10039_c';'R10040_c';'R02985_c';'R02558_c';'R02887_c';'R09375_c';'R09376_c';'R01177_c';'R01896_c';'R08240_c';'R04951_c';'R02569_c';'R00432_c';'R00727_c';'R03664_c';'R11896_c';'R02739_c';'R07346_c';'R02933_c';'R03751_c';'R00802_c';'R06088_c';'R11217_c';'R11219_c';'R04930_c';'R04770_c';'R09366_c';'R00782_c';'R02408_c';'R08193_c';'R00256_c';'R01395_c';'R00575_c';'R10948_c';'R10949_c';'R00717_c';'R01388_c';'R05922_c';'R08549_c';'R00548_c';'R05919_c';'R01000_c';'R03104_c';'R10563_c';'R10147_c';'R01512_c';'R07324_c';'R01993_c';'R03546_c';'R10079_c';'R09365_c';'R05086_c';'R10124_c';'R01078_c';'R01718_c';'R00801_c';'r726';'R06199_c';'R00985_c';'R00236_c';'R00316_c';'R00925_c';'R00926_c';'R01354_c';'R04880_c';'R05233_c';'R05234_c';'R06917_c';'R06927_c';'R07105_c';'R08281_c';'R08306_c';'R08310_c';'R02124_c';'R00754_c';'R00094_c';'R01776_c';'R01777_c';'R00425_c';'r503';'r504';'r505';'r506';'r1398';'r2083';'r564';'r353';'r1743';'r1903';'r1910';'r1917';'r1924';'r307';'r1518';'r1535';'r687';'r930';'r934';'r1941';'r680';'r676';'r1491';'r1667';'r996';'r855';'r193';'r351';'r718';'r570';'r573';'r576';'r579';'r868';'r527';'r1454';'r206';'r660';'r806';'r788';'r793';'r376';'r538';'r1664';'r1451';'r781';'r782';'r783';'r595';'r571';'r574';'r577';'r580';'r728';'r592';'r133';'r873';'r216';'r1523';'r812';'r817';'r2095';'r103';'r1780';'r623';'r112';'r771';'r91';'r878';'r889';'r1938';'r2161';'r90';'r989';'r140';'r141';'r501';'r394';'r1713';'r1831';'r1615';'r1541';'r1653';'r407';'r1860';'r533';'r1796';'R04360_c';'R00286_c';'R07676_c';'R06622_c';'R03566_c';'R06519_c';'R02204_c';'R11165_c';'R11220_c';'R01169_c';'R01989_c';'R04811_c';'R07769_c';'R09419_c';'R10587_c';'R06527_c';'R01461_c';'R06128_c';'R09676_c';'R06238_c';'R06261_c';'R10209_c';'R05976_c';'R07273_c';'R10951_c';'R01059_c';'R03819_c';'R10548_c';'R01547_c';'R04247_c';'R02590_c';'R01618_c';'R09726_c';'R03105_c';'R07461_c';'R07768_c';'R02687_c';'R03416_c';'R03417_c';'R11680_c';'R09827_c';'R03346_c';'R10952_c';'R02541_c';'R01470_c';'R06722_c';'R05981_c';'R06200_c';'R06516_c';'R01647_c';'R00470_c';'R10092_c';'R10827_c';'R01055_c';'R04773_c';'R01724_c';'R03905_c';'R08107_c';'R09073_c';'R10633_c';'r1568';'r978';'r981';'r1313';'r1328';'r2012';'r2064';'R10120_c';'R07763_c';'R10788_c';'R10995_c';'R10851_c';'R10852_c';'R09449_c';'R10828_c';'R03643_c';'R03816_c';'R04866_c';'R04867_c';'R07620_c';'R11399_c';'R01926_c';'R02052_c';'R07381_c';'R02051_c';'R07385_c';'R03980_c';'R04856_c';'R01421_c';'R01515_c';'R04804_c';'R07484_c';'R07771_c';'R11143_c';'R10822_c';'R10823_c';'r1844';'r1849';'r1854';'r1602';'r1616';'r1629';'r1983';'r1992';'r2018';'R07486_c';'R07491_c';'R07505_c';'R01456_c';'R07487_c';'R07492_c';'R02497_c';'R08954_c';'R10242_c';'R02661_c';'R00924_c';'R04432_c';'R02613_c';'R02382_c';'R04300_c';'R00277_c';'R00278_c';'R01710_c';'R04920_c';'R06364_c';'R07384_c';'R09413_c';'R09414_c';'R09415_c';'R06915_c';'R06936_c';'R06939_c';'R05632_c';'R00731_c';'R02363_c';'R04693_c';'R04884_c';'R05800_c';'R05801_c';'R10065_c';'R10953_c';'R11892_c';'R00512_c';'R01665_c';'R00158_c';'R04534_c';'R04536_c';'R04543_c';'R04566_c';'R04953_c';'R04258_c';'R05299_c';'R08114_c';'R08115_c';'R09134_c';'R09035_c';'R09036_c';'R01318_c';'R04480_c';'R03437_c';'R02352_c';'R02353_c';'R04681_c';'R04682_c';'R08945_c';'R08980_c';'R01278_c';'R03776_c';'R03856_c';'R03989_c';'R04753_c';'R06985_c';'R01175_c';'R01279_c';'R03990_c';'R03777_c';'R03857_c';'R04754_c';'R01317_c';'R02053_c';'R07064_c';'R07379_c';'R07387_c';'R07859_c';'r1645';'r1646';'r1647';'r1648';'r1649';'r1650';'R02920_c';'R03304_c';'R04301_c';'R04762_c';'R04764_c';'R04881_c';'R04887_c';'R05333_c';'R07493_c';'R11096_c';'R01457_c';'R03689_c';'R05703_c';'R07498_c';'R07499_c';'R07507_c';'r2003';'r2023';'r2027';'r2031';'r2040';'r2070';'r2052';'r2056';'r2060';'R05510_c';'R05511_c';'R06835_c';'R06838_c';'R08120_c';'R08121_c';'R09136_c';'R09220_c';'R09222_c';'R01103_c';'R01104_c';'R01194_c';'R01329_c';'R02926_c';'R03634_c';'R04019_c';'R04470_c';'R05549_c';'R06091_c';'R05961_c';'R02908_c';'R02919_c';'R04025_c';'R04890_c';'R04893_c';'R04894_c';'R04907_c';'R08346_c';'R08347_c';'R08348_c';'R02173_c';'R04674_c';'R04137_c';'R04204_c';'R04224_c';'R05595_c';'R06942_c';'R03045_c';'R02685_c';'R04170_c';'R04738_c';'R04740_c';'R04744_c';'R04746_c';'R04749_c';'R07002_c';'R07003_c';'R07004_c';'R07023_c';'R07024_c';'R07025_c';'R07026_c';'R07069_c';'R07070_c';'R07083_c';'R07084_c';'R07091_c';'R07092_c';'R07093_c';'R07094_c';'R07100_c';'R07113_c';'R07116_c';'R08280_c';'R09409_c';'R11905_c'};
reducedModel=removeReactions(refineModel,duplicateRxns,true,...
          true,true);
% unbound some reactions
irreversibleRxnsunbound = {'R00078_c';'R00100_c';'R00102_c';'R00126_c';'R00127_c';'R00132_c';'R00196_c';'R00282_c';'R00317_c';'R00381_c';'R00382_c';'R00471_c';'R00565_c';'R00597_c';'R00608_c';'R00610_c';'R00815_c';'R00817_c';'R00818_c';'R00821_c';'R00866_c';'R00966_c';'R01005_c';'R01025_c';'R01030_c';'R01056_c';'R01274_c';'R01280_c';'R01295_c';'R01312_c';'R01315_c';'R01352_c';'R01408_c';'R01416_c';'R01451_c';'R01462_c';'R01494_c';'R01514_c';'R01560_c';'R01626_c';'R01645_c';'R01708_c';'R01711_c';'R01827_c';'R01931_c';'R02025_c';'R02027_c';'R02057_c';'R02144_c';'R02208_c';'R02231_c';'R02239_c';'R02241_c';'R02250_c';'R02323_c';'R02383_c';'R02422_c';'R02480_c';'R02529_c';'R02532_c';'R02534_c';'R02556_c';'R02656_c';'R02746_c';'R02747_c';'R02756_c';'R02964_c';'R03026_c';'R03057_c';'R03172_c';'R03308_c';'R03315_c';'R03332_c';'R03353_c';'R03409_c';'R03435_c';'R03451_c';'R03478_c';'R03522_c';'R03659_c';'R03719_c';'R03772_c';'R03789_c';'R03805_c';'R03828_c';'R03875_c';'R03893_c';'R03938_c';'R03942_c';'R04072_c';'R04095_c';'R04212_c';'R04216_c';'R04404_c';'R04420_c';'R04452_c';'R04496_c';'R04533_c';'R04591_c';'R04620_c';'R04734_c';'R04751_c';'R04756_c';'R04858_c';'R04929_c';'R04964_c';'R05202_c';'R05287_c';'R05556_c';'R05616_c';'R05802_c';'R05916_c';'R05917_c';'R05918_c';'R05920_c';'R05923_c';'R05924_c';'R05970_c';'R05972_c';'R05973_c';'R05979_c';'R05980_c';'R05982_c';'R05986_c';'R06127_c';'R06171_c';'R06258_c';'R06259_c';'R06260_c';'R06262_c';'R06263_c';'R06264_c';'R06513_c';'R06517_c';'R06518_c';'R06601_c';'R07162_c';'R07215_c';'R07459_c';'R07488_c';'R07495_c';'R07758_c';'R07759_c';'R07760_c';'R07761_c';'R07762_c';'R07766_c';'R07767_c';'R07770_c';'R07816_c';'R08266_c';'R08555_c';'R08599_c';'R08701_c';'R09030_c';'R09034_c';'R09037_c';'R09038_c';'R09072_c';'R09087_c';'R09106_c';'R09289_c';'R09304_c';'R09395_c';'R09412_c';'R09597_c';'R09735_c';'R09782_c';'R09844_c';'R09845_c';'R09944_c';'R10116_c';'R10151_c';'R10231_c';'R10309_c';'R10455_c';'R10463_c';'R10532_c';'R10550_c';'R10565_c';'R10586_c';'R10648_c';'R10685_c';'R10686_c';'R10722_c';'R10815_c';'R10826_c';'R11014_c';'R11023_c';'R11166_c';'R11180_c';'R11216_c';'R11218_c';'R11221_c';'R11308_c';'R11583_c';'R11671_c';'R11861_c';'R11891_c';'R11929_c';'R12085_c';'r1113';'r1283';'r1487';'r1495';'r1510';'r1512';'r1526';'r1528';'r1536';'r1546';'r1559';'r1560';'r1561';'r1563';'r1564';'r1569';'r1570';'r1572';'r1576';'r1577';'r1584';'r1587';'r1588';'r1589';'r1644';'r1651';'r171';'r1718';'r1719';'r1727';'r1751';'r1756';'r1784';'r1786';'r1797';'r1800';'r1808';'r1809';'r1814';'r1827';'r1829';'r1839';'r1870';'r1877';'r1988';'r199';'r1997';'r2007';'r203';'r2035';'r204';'r2044';'r2074';'r2093';'r2099';'r2101';'r269';'r388';'r395';'r419';'r427';'r440';'r441';'r443';'r447';'r466';'r473';'r477';'r48';'r502';'r507';'r508';'r509';'r532';'r535';'r561';'r64';'r734';'r765';'r777';'r839';'r867';'r923';'r925';'r975';'r730';'R00533_c';'r1399';'r697';'r668';'R01009_c';'R00086_c';'r499'}; 
unboundModel=setParam(reducedModel,'ub',irreversibleRxnsunbound,1000);
unboundModel=setParam(unboundModel,'lb',irreversibleRxnsunbound,0);

%% STEP 5: MODEL VALIDATION
% The model validation was performed against experimentation.
validateModel = unboundModel;
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
validateModel.description = 'C. militaris GEM by Nachon Raethong [02-AUG-2019]';
validateModel.annotation = [];
exportForGit(validateModel,'model','C:\Users\Dell\Documents\GitHub\Cordyceps_militaris-GEM\');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
