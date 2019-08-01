%To expand curatedModel with iWV1346 and KEGG models

cd 'C:\Users\Dell\Documents\GitHub\cordyceps_model';
aor = importModel('ComplementaryData\iWV1346.xml');
ccm = importModel('ModelFiles\xml\modelATP.xml');
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

%merged = mergeModels({ccm,cmtDraftFromAoryzae});
%exportToExcelFormat(new,'ccmAornew.xlsx');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%cmtKEGG = getKEGGModelForOrganism('cmt');
%save ('ComplementaryData\cmtKEGG.mat','cmtKEGG')
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelBigger = contractModel(newSubSystemModel);
sol = solveLP(modelBigger,1);
fprintf(['Yield of biomass (umax) is '  num2str(sol.f*-1) ' /h\n']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelBigger.id = 'iNR1320';
modelBigger.description = 'Genome-scale metabolic model of C militaris by Nachon [18-MAY-2019]';
modelBigger.annotation = [];
exportForGit(modelBigger,'modelBigger','C:\Users\Dell\Documents\GitHub\cordyceps_model\');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all