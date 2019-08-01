%Reconstruction process
% Purposed: This is the first script used for reconstructing a metabolic model of Cordyceps militaris from yeast model.
% 
% Written by Nachon Raethong, 14-DEC-2018
% Re-upload, 31-JAN-2019
% 
%% WORKSPACE
cd 'C:\Users\Nachonase\Documents\GitHub\Cordyceps_militaris-GEM';
%% stage 1.1:call_addRxns
%% List reactions to add
% 1. biomass 7 rxns
cmtBiomass.rxns = {'cmt_protein';'cmt_DNA';'cmt_RNA';'cmt_lipid';'cmt_carbohydrate';'cmt_vitamins';'cmt_biomass'};
cmtBiomass.rxnNames = {'cmtProtein pseudoreaction';'cmtDNA pseudoreaction';'cmtRNA pseudoreaction';'cmtlipid pseudoreaction';'cmtCarbohydrate pseudoreaction';'cmtVitamins pseudoreaction';'cmtBiomass pseudoreaction'};
cmtBiomass.equations = {'0.28935 L-alanine + 0.18756 L-arginine + 0.09830 L-asparagine + 0.16889 L-aspartate + 0.03792 L-cysteine + 0.16510 L-glutamate + 0.11668 L-glutamine + 0.20447 L-glycine + 0.07263 L-histidine + 0.13009 L-isoleucine + 0.26281 L-leucine + 0.13009 L-lysine + 0.06475 L-methionine + 0.10559 L-phenylalanine + 0.17793 L-proline + 0.22927 L-serine + 0.17297 L-threonine + 0.04259 L-tryptophan + 0.07701 L-tyrosine + 0.18318 L-valine => cmtProtein';...
'0.00507 dAMP + 0.00507 dTMP + 0.00536 dCMP + 0.00536 dGMP => cmtDNA';...
'0.03577 AMP + 0.03577 UMP + 0.04549 CMP + 0.04549 GMP => cmtRNA';...
'0.00016 laurate + 0.00016 myristate + 0.00055 Pentadecanoic acid + 0.07058 palmitate + 0.00564 palmitoleate + 0.02478 oleate + 0.07250 stearate + 0.09821 Linoleate + 0.00197 (9Z,12Z,15Z)-Octadecatrienoic acid + 0.00037 arachidate + 0.00021 (4Z,7Z,10Z,13Z,16Z,19Z)-Docosahexaenoic acid + 0.00026 lignoceric acid + 0.00003 phytosphingosine + 0.00002 sphinganine + 0.00001 sphinganine 1-phosphate + 0.00003 Sphingomyelin + 0.00004 Phytoceramide + 0.00078 ergosterol + 0.05253 glycerol 3-phosphate + 0.01086 CDP-choline + 0.00757 CDP-ethanolamine => cmtLipid';...
'0.63086 (1->3)-beta-D-glucan + 0.02465 UDP-N-acetyl-alpha-D-glucosamine + 0.02024 GDP-alpha-D-mannose + 0.11020 UDP-D-glucose + 0.01437 UDP-D-galactose => cmtCarbohydrate';...
'0.00809 Retinol + 0.01321 Ergothioneine + 0.16830 Cordycepin + 0.000000001 Deoxycoformycin + 0.00662 Mannitol + 0.04049 gamma-aminobutyrate + 0.02 sulphate + 0.00002 riboflavin => cmtVitamins';...
'88.0022 ATP + 88.0022 H2O + cmtProtein + cmtDNA + cmtRNA + cmtCarbohydrate + cmtLipid + cmtVitamins => cmtBiomass + 88.0022 ADP + 88.0022 phosphate + 88.0022 H+'};


% 2. sphingo lipid 14 rxns
cmtSphingo.rxns = {'cmt_ceramide-1';'cmt_ceramide-2';'cmt_ceramide-3';'cmt_ceramide-4';'cmt_ceramide-5';'cmt_ceramide-6';'cmt_ceramide-7';'cmt_ceramide-8';'cmt_ceramide-9';'cmt_ceramide-10';'cmt_sphingomyelin-1';'cmt_ethanolamine-1';'cmt_ethanolamine-2';'cmt_ethanolamine-3'};
cmtSphingo.equations = {'sphinganine + tetracosanoyl-CoA => Phytoceramide + coenzyme A + H+';...
'hexacosanoyl-CoA + sphinganine => Phytoceramide + coenzyme A + H+';...
'phytosphingosine + tetracosanoyl-CoA => Phytoceramide + coenzyme A + H+';...
'hexacosanoyl-CoA + phytosphingosine => Phytoceramide + coenzyme A + H+';...
'lignoceric acid + sphinganine => Phytoceramide';...
'cerotic acid + sphinganine => Phytoceramide';...
'Phytoceramide => lignoceric acid + phytosphingosine';...
'Phytoceramide => cerotic acid + phytosphingosine';...
'lignoceric acid + phytosphingosine => Phytoceramide';...
'cerotic acid + phytosphingosine => Phytoceramide';...
'Phytoceramide + CDP-choline => Sphingomyelin + CMP + H+';...
'phosphate + acetaldehyde + ammonium => H2O + O-phosphoethanolamine';...
'hydrogen peroxide + glycolaldehyde + ammonium => oxygen + H2O + ethanolamine';...
'L-serine => carbon dioxide + ethanolamine'};

% 3. long chain lipid 7 rxns
cmtLongchain.rxns = {'cmt_docosahexaenoate-1';'cmt_linoleate-1';...
    'cmt_linoleate-2';'cmt_linolenate-1';'cmt_linolenate-2';...
    'cmt_pentadecanoate-1';'cmt_pentadecanoate-2'};
cmtLongchain.equations = {'acetyl-CoA + 23 H+ + 10 malonyl-CoA + 14 NADPH => 10 carbon dioxide + 11 coenzyme A + 9 H2O + 14 NADP(+) + (4Z,7Z,10Z,13Z,16Z,19Z)-Docosahexaenoic acid';...
    'oleoyl-CoA + oxygen + 2 H+ => Linoleoyl-CoA + + 2 H2O';...
    'Linoleoyl-CoA + H2O => coenzyme A + Linoleate';...
    'Linoleoyl-CoA + oxygen + 2 H+ => (9Z,12Z,15Z)-Octadecatrienoyl-CoA + 2 H2O';...
    '(9Z,12Z,15Z)-Octadecatrienoyl-CoA + H2O => coenzyme A + (9Z,12Z,15Z)-Octadecatrienoic acid';...
    'palmitate + 2 hydrogen peroxide => Pentadecanal + carbon dioxide + 3 H2O';...
    'Pentadecanal + NAD + H2O => Pentadecanoic acid + NADH + 2 H+'};

% 4. mannitol 3 rxns
cmtMannitol.rxns = {'cmt_mannitol-1';'cmt_mannitol-2';'cmt_mannitol-3'};
cmtMannitol.equations = {'NADH + D-fructose 6-phosphate => NAD + D-Mannitol 1-phosphate + H+';...
    'D-mannose 6-phosphate + NADPH + H+ => D-Mannitol 1-phosphate + NADP(+)';...
    'D-Mannitol 1-phosphate + H2O => Mannitol + diphosphate'};
	

% 5. cordycepin 6 rxns
cmtCordycepin.rxns = {'cmt_cordycepin-1';'cmt_cordycepin-2';'cmt_cordycepin-3';...
    'cmt_cordycepin-4';'cmt_cordycepin-5';'cmt_cordycepin-6'};
cmtCordycepin.equations = {'ATP + cmtRNA <=> 2'',3''-Cyclic AMP + diphosphate + RNA(circular)';...
    '2'',3''-Cyclic AMP => 3''-AMP + Adenosine 2''-phosphate';...
    'adenosine + ATP => ADP + phosphate + 3''-AMP';...
    '3''-AMP => Cordycepin';...
    'Cordycepin => 2''-deoxyadenosine';...
    'adenosine + ATP + PRPP => Deoxycoformycin'};

% 6. ergothioneine 3 rxns
cmtErgothioneine.rxns = {'cmt_ergothioneine-1';'cmt_ergothioneine-2';'cmt_ergothioneine-3'};
cmtErgothioneine.equations = {'L-histidine + 3 S-adenosyl-L-methionine + H+ => Hercynine + 3 S-adenosyl-L-homocysteine';...
    'Hercynine + L-cysteine + oxygen => S-(Hercyn-2-yl)-L-cysteine S-oxide + H2O';...
    'S-(Hercyn-2-yl)-L-cysteine S-oxide => Ergothioneine + pyruvate'};


% 7. retinol 12 rxns
cmtRetinol.rxns = {'cmt_retinol-1';'cmt_retinol-2';'cmt_retinol-3';'cmt_retinol-4';...
    'cmt_retinol-5';'cmt_retinol-6';'cmt_retinol-7';'cmt_retinol-8';...
    'cmt_retinol-9';'cmt_retinol-10';'cmt_retinol-11';'cmt_retinol-12'};
cmtRetinol.equations = {'geranylgeranyl diphosphate => Prephytoene diphosphate + diphosphate';...
    'Prephytoene diphosphate => 15-cis-Phytoene + diphosphate';...
    '15-cis-Phytoene => 15,9''-dicis-Phytofluene';...
    '15,9''-dicis-Phytofluene => 9,15,9''-tricis-zeta-Carotene';...
    '9,15,9''-tricis-zeta-Carotene => 9,9''-dicis-zeta-Carotene';...
    '9,9''-dicis-zeta-Carotene => 7,9,9''-tricis-Neurosporene';...
    '7,9,9''-tricis-Neurosporene => 7,9,7'',9''-tetracis-Lycopene';...
    '7,9,7'',9''-tetracis-Lycopene => Lycopene';...
    'Lycopene => gamma-Carotene';...
    'gamma-Carotene => beta-Carotene';...
    'beta-Carotene => Retinal';...
    'NADH + H+ + Retinal <=> Retinol + NAD'};

% 8. carbon source utilization 5 rxns
cmtUtilization.rxns = {'cmt_lactose-1';'cmt_lactose-2';...
    'cmt_arabinose-1';'cmt_arabinose-2';...
    'cmt_cellobiose-1'};
cmtUtilization.equations = {'Lactose + H2O => D-glucose + D-galactose';...
    'UDP + Lactose => UDP-D-galactose + D-glucose';...
    'ATP + D-arabinose => D-Arabinose 5-phosphate + ADP';...
    'D-Arabinose 5-phosphate <=> D-ribulose 5-phosphate';...
    'Cellobiose + diphosphate <=> D-glucose 1-phosphate + D-glucose'};

% 9. exchange rxns 19 rxns (_ex_c_:11 + _ex_:8)
cmtExchange.rxns = {'cmt_ex_c_glucose';...
    'cmt_ex_c_arabinose';...
    'cmt_ex_c_fructose';...
    'cmt_ex_c_sorbose';...
    'cmt_ex_c_mannose';...
    'cmt_ex_c_xylose';...
    'cmt_ex_c_sucrose';...
    'cmt_ex_c_trehalose';...
    'cmt_ex_c_maltose';...
    'cmt_ex_c_lactose';...
    'cmt_ex_c_cellobiose'};
cmtExchange.equations = {'D-glucose <=> ';...
    'D-arabinose <=> ';...
    'D-fructose <=> ';...
    'D-glucitol <=> ';...
    'D-mannose <=> ';...
    'D-xylose <=> ';...
    'sucrose <=> ';...
    'trehalose <=> ';...
    'maltose <=> ';...
    'Lactose <=> ';...
    'Cellobiose <=> '};

% 10. objective function 2 rxns
cmtObjective.rxns = {'cmt_ex_growth';'cmt_ex_cordycepin'};
cmtObjective.rxnNames = {'C. militaris growth';'extracellular cordycepin'};
cmtObjective.equations = {'cmtBiomass => ';'Cordycepin => '};

% pool 78 rxns
cmt_addRxns.rxns = [cmtBiomass.rxns;...
    cmtSphingo.rxns;...
    cmtLongchain.rxns;...
    cmtMannitol.rxns;...
    cmtCordycepin.rxns;...
    cmtRetinol.rxns;...
    cmtErgothioneine.rxns;...
    cmtUtilization.rxns;...
    cmtExchange.rxns;...
    cmtObjective.rxns];

cmt_addRxns.equations = [cmtBiomass.equations;...
    cmtSphingo.equations;...
    cmtLongchain.equations;...
    cmtMannitol.equations;...
    cmtCordycepin.equations;...
    cmtRetinol.equations;...
    cmtErgothioneine.equations;...
    cmtUtilization.equations;...
    cmtExchange.equations;...
    cmtObjective.equations];

%% stage 1.2: call_tmpModel
load 'ComplementaryData\yeastGEM_v8.0.2.mat\'; % obtained from SysBioChalmers github (25-MAY-2018)

yeastGEM = ravenCobraWrapper(model); % convert yeastGEM from cobra to raven format ..
yeastGEM.description = 'yeastGEM_v8.0.2';
yeastGEM.equations = constructEquations(yeastGEM);
tmpModel = mergeCompartments(yeastGEM); % yeast model with 1 compartment
tmpModel.equations = constructEquations(tmpModel);
tmpModel.id = 'yeast'; % must same with blastStructure {'yeast'}
tmpModel.description = 'template model';
tmpModel.rxnKEGGID = cell(numel(tmpModel.rxns),1);

for i = 1:numel(tmpModel.rxns)
    I = find(ismember(model.rxns,tmpModel.rxns(i)));
    tmpModel.rxnKEGGID(i) = model.rxnKEGGID(I);
end

tmpModel=setParam(tmpModel,'eq','r_1714',-1); % glucose uptake
tmpModel=setParam(tmpModel,'obj','r_2111',1);
solveLP(tmpModel);
sol = solveLP(tmpModel);
tmpModel.sol = sol.x;
printFluxes(tmpModel,sol.x,true);

%% stage 2: call_cmtCore
load('C:\Users\Nachonase\Documents\GitHub\Cordyceps_militaris-GEM\ComplementaryData\keggMets.mat');
[~,metNames] = xlsread('C:\Users\Nachonase\Documents\GitHub\Cordyceps_militaris-GEM\ComplementaryData\metsToAdd.xlsx','metNames'); 
[~,mets] = xlsread('C:\Users\Nachonase\Documents\GitHub\Cordyceps_militaris-GEM\ComplementaryData\metsToAdd.xlsx','mets'); 
[~,metFormulas] = xlsread('C:\Users\Nachonase\Documents\GitHub\Cordyceps_militaris-GEM\ComplementaryData\metsToAdd.xlsx','metFormulas'); 

metsToAdd = struct();
metsToAdd.mets = mets(1:end,1);
metsToAdd.metNames = metNames(1:end,1);
metsToAdd.metFormulas = metFormulas(1:end,1);
metsToAdd.inchis = cell(numel(metsToAdd.metNames),1);
metsToAdd.metMiriams = cell(numel(metsToAdd.metNames),1);
metsToAdd.compartments = cell(numel(metsToAdd.metNames),1);
for i=1:numel(metsToAdd.metNames)
    metsToAdd.compartments{i} = 's';
end
[a, b] = ismember(metsToAdd.metNames,keggMets.metNames);
I = find(a);
metsToAdd.inchis(I) = keggMets.inchis(b(I));
metsToAdd.metMiriams(I) = keggMets.metMiriams(b(I));

dummyModel=addMets(tmpModel,metsToAdd);


cmtRxns = addRxns(dummyModel,cmt_addRxns,2,'s',true);
cmtRxns.equations = constructEquations(cmtRxns);

cmtRxns = removeReactions(cmtRxns,{'r_1650';'r_1706';'r_1709';'r_1712';...
    'r_1715';'r_1718';'r_1931';'r_2058';...
    'r_1714'},true);


id_cmtRxns = regexp(cmtRxns.rxns,...
    'cmt_','all');
cmt_addRxns = cmtRxns.rxns(find(~cellfun('isempty',id_cmtRxns)));

for i = 1:numel(cmt_addRxns)
    if cmtRxns.rev(find(ismember(cmtRxns.rxns,cmt_addRxns(i)))) == 0;
        cmtRxns.lb(find(ismember(cmtRxns.rxns,cmt_addRxns(i)))) = 0;
        cmtRxns.ub(find(ismember(cmtRxns.rxns,cmt_addRxns(i)))) = 1000;
    else 
        cmtRxns.lb(find(ismember(cmtRxns.rxns,cmt_addRxns(i)))) = -1000;
        cmtRxns.ub(find(ismember(cmtRxns.rxns,cmt_addRxns(i)))) = 1000;
    end
end
    
%ATP leak rxns
cmtRxns=setParam(cmtRxns,'eq',{'r_1022'},0);
cmtRxns=setParam(cmtRxns,'eq',{'r_2212'},0);

cmtRxns=setParam(cmtRxns,'lb','cmt_ex_c_glucose',-1);
cmtRxns=setParam(cmtRxns,'obj','cmt_ex_growth',1);
sol = solveLP(cmtRxns);
printFluxes(cmtRxns,sol.x,true);

%%
id_cmt_ex_c = regexp(cmtRxns.rxns,...
    'cmt_ex_c_','all');
cmt_ex_c = cmtRxns.rxns(find(~cellfun('isempty',id_cmt_ex_c)));
cmtRxns=setParam(cmtRxns,'lb','cmt_ex_cordycepin',0.0000001);
cmtRxns=setParam(cmtRxns,'obj','cmt_ex_growth',1);
sol = solveLP(cmtRxns);
cmtRxns.haveFlux(1) = sol;
for i = 1:numel(cmt_ex_c)
    cmtRxns=setParam(cmtRxns,'lb',cmt_ex_c,0);
    cmtRxns=setParam(cmtRxns,'lb',cmt_ex_c(i),-1);
    sol = solveLP(cmtRxns);
    cmtRxns.haveFlux(i) = sol;
end

cmtRxns.compareFlux = table(cmtRxns.rxns,cmtRxns.equations,cmtRxns.haveFlux(1).x,...
    cmtRxns.haveFlux(2).x,cmtRxns.haveFlux(3).x,cmtRxns.haveFlux(4).x,...
    cmtRxns.haveFlux(5).x,cmtRxns.haveFlux(6).x,cmtRxns.haveFlux(7).x,...
    cmtRxns.haveFlux(8).x,cmtRxns.haveFlux(9).x,cmtRxns.haveFlux(10).x,...
    cmtRxns.haveFlux(11).x);

non_essential={'r_0003';'r_0004';'r_0006';'r_0013';'r_0017';'r_0019';'r_0021';'r_0022';'r_0026';'r_0028';'r_0033';'r_0034';'r_0035';'r_0036';'r_0037';'r_0043';'r_0044';'r_0045';'r_0058';'r_0059';'r_0062';'r_0063';'r_0064';'r_0066';'r_0067';'r_0068';'r_0069';'r_0070';'r_0072';'r_0073';'r_0074';'r_0075';'r_0076';'r_0077';'r_0078';'r_0081';'r_0082';'r_0083';'r_0084';'r_0085';'r_0086';'r_0087';'r_0088';'r_0089';'r_0090';'r_0091';'r_0092';'r_0093';'r_0094';'r_0095';'r_0099';'r_0101';'r_0106';'r_0107';'r_0111';'r_0116';'r_0117';'r_0119';'r_0120';'r_0121';'r_0126';'r_0127';'r_0128';'r_0129';'r_0130';'r_0131';'r_0132';'r_0133';'r_0134';'r_0135';'r_0137';'r_0138';'r_0139';'r_0140';'r_0142';'r_0143';'r_0145';'r_0146';'r_0147';'r_0155';'r_0156';'r_0157';'r_0158';'r_0159';'r_0160';'r_0161';'r_0162';'r_0163';'r_0164';'r_0165';'r_0166';'r_0168';'r_0169';'r_0171';'r_0172';'r_0174';'r_0176';'r_0177';'r_0179';'r_0181';'r_0182';'r_0184';'r_0185';'r_0186';'r_0188';'r_0189';'r_0190';'r_0191';'r_0192';'r_0195';'r_0199';'r_0200';'r_0201';'r_0204';'r_0205';'r_0206';'r_0209';'r_0212';'r_0217';'r_0220';'r_0222';'r_0223';'r_0224';'r_0227';'r_0228';'r_0229';'r_0230';'r_0233';'r_0242';'r_0243';'r_0249';'r_0252';'r_0254';'r_0259';'r_0260';'r_0261';'r_0262';'r_0263';'r_0264';'r_0265';'r_0266';'r_0267';'r_0268';'r_0269';'r_0270';'r_0271';'r_0272';'r_0282';'r_0283';'r_0284';'r_0285';'r_0286';'r_0287';'r_0288';'r_0289';'r_0290';'r_0291';'r_0292';'r_0293';'r_0294';'r_0295';'r_0296';'r_0297';'r_0298';'r_0299';'r_0304';'r_0308';'r_0311';'r_0313';'r_0314';'r_0315';'r_0318';'r_0319';'r_0320';'r_0321';'r_0322';'r_0327';'r_0328';'r_0329';'r_0331';'r_0332';'r_0334';'r_0335';'r_0340';'r_0341';'r_0342';'r_0343';'r_0346';'r_0347';'r_0348';'r_0350';'r_0351';'r_0356';'r_0357';'r_0358';'r_0359';'r_0360';'r_0361';'r_0362';'r_0363';'r_0365';'r_0368';'r_0369';'r_0370';'r_0399';'r_0410';'r_0436';'r_0437';'r_0440';'r_0441';'r_0442';'r_0443';'r_0445';'r_0448';'r_0449';'r_0454';'r_0457';'r_0459';'r_0460';'r_0463';'r_0465';'r_0466';'r_0468';'r_0470';'r_0473';'r_0475';'r_0478';'r_0479';'r_0481';'r_0483';'r_0485';'r_0488';'r_0490';'r_0497';'r_0500';'r_0501';'r_0510';'r_0511';'r_0512';'r_0518';'r_0519';'r_0520';'r_0521';'r_0522';'r_0523';'r_0524';'r_0526';'r_0527';'r_0530';'r_0531';'r_0532';'r_0539';'r_0541';'r_0544';'r_0547';'r_0550';'r_0555';'r_0556';'r_0557';'r_0561';'r_0562';'r_0567';'r_0571';'r_0572';'r_0573';'r_0574';'r_0575';'r_0596';'r_0597';'r_0598';'r_0599';'r_0600';'r_0601';'r_0602';'r_0603';'r_0604';'r_0605';'r_0606';'r_0607';'r_0608';'r_0609';'r_0610';'r_0611';'r_0612';'r_0613';'r_0614';'r_0615';'r_0616';'r_0617';'r_0618';'r_0619';'r_0620';'r_0621';'r_0622';'r_0623';'r_0624';'r_0625';'r_0646';'r_0647';'r_0648';'r_0649';'r_0650';'r_0651';'r_0652';'r_0653';'r_0654';'r_0655';'r_0656';'r_0657';'r_0658';'r_0661';'r_0662';'r_0665';'r_0668';'r_0670';'r_0671';'r_0672';'r_0673';'r_0675';'r_0676';'r_0679';'r_0681';'r_0687';'r_0688';'r_0690';'r_0691';'r_0694';'r_0695';'r_0696';'r_0701';'r_0703';'r_0705';'r_0707';'r_0711';'r_0716';'r_0718';'r_0721';'r_0728';'r_0729';'r_0731';'r_0734';'r_0737';'r_0747';'r_0748';'r_0749';'r_0750';'r_0751';'r_0752';'r_0753';'r_0754';'r_0755';'r_0756';'r_0757';'r_0758';'r_0761';'r_0762';'r_0763';'r_0764';'r_0765';'r_0767';'r_0770';'r_0771';'r_0774';'r_0781';'r_0782';'r_0783';'r_0785';'r_0786';'r_0788';'r_0789';'r_0790';'r_0793';'r_0797';'r_0798';'r_0799';'r_0801';'r_0802';'r_0803';'r_0804';'r_0805';'r_0806';'r_0807';'r_0812';'r_0815';'r_0817';'r_0819';'r_0831';'r_0832';'r_0841';'r_0842';'r_0843';'r_0844';'r_0845';'r_0847';'r_0848';'r_0849';'r_0850';'r_0852';'r_0854';'r_0884';'r_0889';'r_0890';'r_0903';'r_0905';'r_0906';'r_0919';'r_0920';'r_0921';'r_0922';'r_0929';'r_0935';'r_0936';'r_0937';'r_0940';'r_0941';'r_0942';'r_0943';'r_0949';'r_0951';'r_0953';'r_0954';'r_0955';'r_0956';'r_0959';'r_0960';'r_0963';'r_0965';'r_0971';'r_0972';'r_0983';'r_0985';'r_0986';'r_0987';'r_0995';'r_0998';'r_0999';'r_1000';'r_1001';'r_1002';'r_1003';'r_1004';'r_1005';'r_1006';'r_1008';'r_1011';'r_1022';'r_1023';'r_1029';'r_1030';'r_1031';'r_1032';'r_1033';'r_1034';'r_1035';'r_1036';'r_1042';'r_1046';'r_1051';'r_1056';'r_1057';'r_1066';'r_1068';'r_1074';'r_1075';'r_1076';'r_1077';'r_1078';'r_1079';'r_1081';'r_1082';'r_1083';'r_1089';'r_1091';'r_1095';'r_1619';'r_2029';'r_2112';'r_2113';'r_2114';'r_2116';'r_2118';'r_2126';'r_2142';'r_2143';'r_2144';'r_2145';'r_2146';'r_2147';'r_2148';'r_2149';'r_2150';'r_2151';'r_2152';'r_2153';'r_2154';'r_2155';'r_2160';'r_2161';'r_2162';'r_2167';'r_2169';'r_2174';'r_2175';'r_2176';'r_2210';'r_2212';'r_2232';'r_2233';'r_2234';'r_2235';'r_2236';'r_2237';'r_2238';'r_2239';'r_2240';'r_2242';'r_2243';'r_2244';'r_2245';'r_2246';'r_2247';'r_2248';'r_2249';'r_2250';'r_2252';'r_2253';'r_2254';'r_2255';'r_2256';'r_2257';'r_2258';'r_2259';'r_2260';'r_2261';'r_2262';'r_2263';'r_2264';'r_2265';'r_2266';'r_2267';'r_2270';'r_2271';'r_2272';'r_2273';'r_2274';'r_2275';'r_2276';'r_2277';'r_2278';'r_2279';'r_2280';'r_2281';'r_2282';'r_2284';'r_2285';'r_2286';'r_2287';'r_2288';'r_2289';'r_2290';'r_2291';'r_2292';'r_2293';'r_2294';'r_2295';'r_2296';'r_2297';'r_2298';'r_2299';'r_2300';'r_2301';'r_2302';'r_2303';'r_2304';'r_2308';'r_2309';'r_2310';'r_2311';'r_2312';'r_2313';'r_2314';'r_2315';'r_2324';'r_2325';'r_2326';'r_2327';'r_2332';'r_2333';'r_2334';'r_2335';'r_2336';'r_2337';'r_2338';'r_2339';'r_2344';'r_2349';'r_2368';'r_2369';'r_2370';'r_2371';'r_2372';'r_2373';'r_2374';'r_2375';'r_2376';'r_2377';'r_2380';'r_2381';'r_2382';'r_2383';'r_2384';'r_2385';'r_2386';'r_2387';'r_2388';'r_2389';'r_2390';'r_2391';'r_2392';'r_2393';'r_2394';'r_2395';'r_2396';'r_2397';'r_2399';'r_2432';'r_2433';'r_2435';'r_2436';'r_2437';'r_2438';'r_2439';'r_2446';'r_2447';'r_2448';'r_2449';'r_2450';'r_2451';'r_2452';'r_2453';'r_2454';'r_2455';'r_2456';'r_2457';'r_2458';'r_2459';'r_2460';'r_2461';'r_2462';'r_2463';'r_2464';'r_2465';'r_2466';'r_2467';'r_2468';'r_2469';'r_2470';'r_2471';'r_2489';'r_2490';'r_2491';'r_2494';'r_2497';'r_2498';'r_2499';'r_2502';'r_2505';'r_2506';'r_2507';'r_2510';'r_2512';'r_2515';'r_2516';'r_2518';'r_2519';'r_2520';'r_2522';'r_2526';'r_2528';'r_2529';'r_2530';'r_2531';'r_2532';'r_2533';'r_2534';'r_2535';'r_2536';'r_2537';'r_2539';'r_2540';'r_2541';'r_2542';'r_2543';'r_2545';'r_2546';'r_2547';'r_2548';'r_2549';'r_2550';'r_2551';'r_2552';'r_2553';'r_2554';'r_2555';'r_2556';'r_2557';'r_2558';'r_2559';'r_2560';'r_2561';'r_2562';'r_2563';'r_2564';'r_2565';'r_2566';'r_2567';'r_2568';'r_2569';'r_2570';'r_2571';'r_2572';'r_2573';'r_2574';'r_2575';'r_2576';'r_2577';'r_2578';'r_2579';'r_2580';'r_2581';'r_2582';'r_2583';'r_2584';'r_2585';'r_2586';'r_2587';'r_2588';'r_2589';'r_2590';'r_2591';'r_2592';'r_2593';'r_2594';'r_2595';'r_2596';'r_2597';'r_2598';'r_2599';'r_2600';'r_2601';'r_2602';'r_2603';'r_2604';'r_2605';'r_2606';'r_2607';'r_2608';'r_2609';'r_2610';'r_2611';'r_2612';'r_2613';'r_2614';'r_2615';'r_2616';'r_2617';'r_2618';'r_2619';'r_2820';'r_2821';'r_2822';'r_2823';'r_2824';'r_2825';'r_2826';'r_2827';'r_2852';'r_2853';'r_2854';'r_2855';'r_2856';'r_2857';'r_2858';'r_2859';'r_2860';'r_2861';'r_2862';'r_2863';'r_2864';'r_2865';'r_2866';'r_2867';'r_2876';'r_2877';'r_2878';'r_2879';'r_2880';'r_2881';'r_2882';'r_2883';'r_3022';'r_3024';'r_3026';'r_3027';'r_3028';'r_3029';'r_3030';'r_3031';'r_3032';'r_3033';'r_3046';'r_3047';'r_3048';'r_3049';'r_3051';'r_3052';'r_3053';'r_3054';'r_3055';'r_3056';'r_3057';'r_3058';'r_3059';'r_3060';'r_3061';'r_3062';'r_3063';'r_3064';'r_3065';'r_3066';'r_3067';'r_3068';'r_3069';'r_3070';'r_3071';'r_3072';'r_3073';'r_3074';'r_3075';'r_3076';'r_3077';'r_3078';'r_3079';'r_3080';'r_3081';'r_3082';'r_3083';'r_3084';'r_3085';'r_3086';'r_3087';'r_3088';'r_3089';'r_3098';'r_3099';'r_3101';'r_3102';'r_3103';'r_3104';'r_3106';'r_3110';'r_3112';'r_3113';'r_3114';'r_3115';'r_3116';'r_3117';'r_3118';'r_3119';'r_3144';'r_3145';'r_3146';'r_3147';'r_3148';'r_3149';'r_3150';'r_3151';'r_3176';'r_3177';'r_3178';'r_3179';'r_3180';'r_3181';'r_3182';'r_3183';'r_3192';'r_3193';'r_3194';'r_3195';'r_3196';'r_3197';'r_3198';'r_3199';'r_3224';'r_3225';'r_3226';'r_3227';'r_3228';'r_3229';'r_3230';'r_3231';'r_3240';'r_3241';'r_3242';'r_3243';'r_3252';'r_3253';'r_3254';'r_3255';'r_3256';'r_3257';'r_3258';'r_3259';'r_3260';'r_3261';'r_3264';'r_3265';'r_3266';'r_3267';'r_3268';'r_3269';'r_3270';'r_3271';'r_3272';'r_3273';'r_3274';'r_3276';'r_3277';'r_3278';'r_3279';'r_3280';'r_3281';'r_3282';'r_3283';'r_3284';'r_3285';'r_3286';'r_3287';'r_3288';'r_3289';'r_3290';'r_3291';'r_3292';'r_3293';'r_3294';'r_3295';'r_3296';'r_3297';'r_3298';'r_3299';'r_3300';'r_3301';'r_3302';'r_3303';'r_3308';'r_3309';'r_3310';'r_3311';'r_3312';'r_3313';'r_3314';'r_3315';'r_4039';'r_4042';'r_3348';'r_3349';'r_3350';'r_3351';'r_3352';'r_3353';'r_3354';'r_3355';'r_3356';'r_3357';'r_3358';'r_3359';'r_3360';'r_3361';'r_3362';'r_3363';'r_3364';'r_3365';'r_3366';'r_3367';'r_3368';'r_3369';'r_3370';'r_3371';'r_3372';'r_3373';'r_3374';'r_3375';'r_3376';'r_3377';'r_3378';'r_3379';'r_3380';'r_3381';'r_3382';'r_3383';'r_3384';'r_3385';'r_3386';'r_3387';'r_3388';'r_3389';'r_3390';'r_3391';'r_3392';'r_3393';'r_3394';'r_3395';'r_3396';'r_3397';'r_3398';'r_3399';'r_3400';'r_3401';'r_3402';'r_3403';'r_3404';'r_3405';'r_3406';'r_3407';'r_3408';'r_3409';'r_3410';'r_3411';'r_3412';'r_3413';'r_3414';'r_3415';'r_3416';'r_3417';'r_3418';'r_3419';'r_3420';'r_3421';'r_3422';'r_3423';'r_3424';'r_3425';'r_3426';'r_3427';'r_3428';'r_3429';'r_3430';'r_3431';'r_3432';'r_3433';'r_3434';'r_3435';'r_3436';'r_3437';'r_3438';'r_3439';'r_3440';'r_3441';'r_3442';'r_3443';'r_3444';'r_3445';'r_3446';'r_3447';'r_3448';'r_3449';'r_3450';'r_3451';'r_3452';'r_3453';'r_3454';'r_3455';'r_3456';'r_3457';'r_3458';'r_3459';'r_3460';'r_3461';'r_3462';'r_3463';'r_3464';'r_3465';'r_3466';'r_3467';'r_3468';'r_3469';'r_3470';'r_3471';'r_3472';'r_3473';'r_3474';'r_3475';'r_3476';'r_3477';'r_3478';'r_3479';'r_3480';'r_3481';'r_3482';'r_3483';'r_3484';'r_3485';'r_3486';'r_3487';'r_3488';'r_3489';'r_3490';'r_3491';'r_3492';'r_3493';'r_3494';'r_3495';'r_3496';'r_3497';'r_3498';'r_3499';'r_3500';'r_3501';'r_3502';'r_3503';'r_3504';'r_3505';'r_3506';'r_3507';'r_1357';'r_1358';'r_1359';'r_1363';'r_1364';'r_1365';'r_1366';'r_1367';'r_1368';'r_1369';'r_1370';'r_1371';'r_1449';'r_1450';'r_1451';'r_1452';'r_1453';'r_1454';'r_1455';'r_1456';'r_1457';'r_1458';'r_1479';'r_1480';'r_1481';'r_1482';'r_1483';'r_1484';'r_1485';'r_1486';'r_1487';'r_1488';'r_1509';'r_1510';'r_1511';'r_1512';'r_1513';'r_1514';'r_1515';'r_1516';'r_1517';'r_1518';'r_3963';'r_3964';'r_3965';'r_3966';'r_3967';'r_3968';'r_3969';'r_3970';'r_3971';'r_3972';'r_3974';'r_3975';'r_3976';'r_3977';'r_3978';'r_3979';'r_3980';'r_3981';'r_3982';'r_3983';'r_3984';'r_3985';'r_3986';'r_3988';'r_3989';'r_3990';'r_3991';'r_3992';'r_3993';'r_3994';'r_3995';'r_3997';'r_3998';'r_3999';'r_4000';'r_4001';'r_4002';'r_4003';'r_4004';'r_4006';'r_4007';'r_4008';'r_4009';'r_4010';'r_4011';'r_4012';'r_4013';'r_4014';'r_4015';'r_4016';'r_4017';'r_4018';'r_4019';'r_4020';'r_4021';'r_4022';'r_4023';'r_4024';'r_4025';'r_4026';'r_4027';'r_4028';'r_4029';'r_4030';'r_4031';'r_4032';'r_4033';'r_4034';'r_4035';'r_4036';'r_4037';'r_1542';'r_1545';'r_1546';'r_1547';'r_1548';'r_1549';'r_1550';'r_1551';'r_1552';'r_1553';'r_1554';'r_1563';'r_1564';'r_1565';'r_1566';'r_1572';'r_1577';'r_1580';'r_1581';'r_1586';'r_1589';'r_1598';'r_1603';'r_1604';'r_1614';'r_1615';'r_1616';'r_1617';'r_1618';'r_1620';'r_1621';'r_1624';'r_1625';'r_1627';'r_1629';'r_1630';'r_1631';'r_1634';'r_1639';'r_1641';'r_1643';'r_1648';'r_1649';'r_1651';'r_1671';'r_1683';'r_1685';'r_1687';'r_1690';'r_1702';'r_1703';'r_1705';'r_1710';'r_1711';'r_1716';'r_1722';'r_1723';'r_1724';'r_1725';'r_1727';'r_1730';'r_1734';'r_1739';'r_1744';'r_1749';'r_1753';'r_1757';'r_1761';'r_1764';'r_1765';'r_1788';'r_1791';'r_1792';'r_1793';'r_1797';'r_1798';'r_1800';'r_1806';'r_1807';'r_1808';'r_1810';'r_1814';'r_1815';'r_1818';'r_1820';'r_1821';'r_1834';'r_1841';'r_1843';'r_1847';'r_1850';'r_1861';'r_1862';'r_1865';'r_1866';'r_1867';'r_1870';'r_1871';'r_1872';'r_1873';'r_1875';'r_1878';'r_1879';'r_1880';'r_1881';'r_1883';'r_1886';'r_1887';'r_1889';'r_1891';'r_1893';'r_1896';'r_1897';'r_1899';'r_1900';'r_1902';'r_1903';'r_1904';'r_1906';'r_1909';'r_1911';'r_1912';'r_1913';'r_1914';'r_1915';'r_1916';'r_1947';'r_1952';'r_1967';'r_1968';'r_1974';'r_1975';'r_1984';'r_1987';'r_1989';'r_1993';'r_1994';'r_1999';'r_2000';'r_2001';'r_2020';'r_2024';'r_2025';'r_2026';'r_2027';'r_2028';'r_2033';'r_2038';'r_2043';'r_2044';'r_2046';'r_2049';'r_2050';'r_2051';'r_2052';'r_2055';'r_2056';'r_2061';'r_2062';'r_2065';'r_2066';'r_2067';'r_2068';'r_2069';'r_2070';'r_2071';'r_2073';'r_2074';'r_2083';'r_2090';'r_2091';'r_2092';'r_2102';'r_2104';'r_2106';'r_2108';'r_2111';'r_2134';'r_2137';'r_2187';'r_2188';'r_2189';'r_2193';'r_2812';'r_2813';'r_2814';'r_2815';'r_2816';'r_2817';'r_2818';'r_2819';'r_3332';'r_3333';'r_3334';'r_3335';'r_3336';'r_3337';'r_3338';'r_3339';'r_4043';'r_4044';'r_4041';'r_4047';'r_4048';'r_4049';'r_4050'};
cmtCore = removeReactions(cmtRxns,non_essential,true);
cmtCore=setParam(cmtCore,'obj','cmt_ex_growth',1);
solveLP(cmtCore);

%% stage 3: call_cmtModel
load 'C:\Users\Nachonase\Documents\GitHub\Cordyceps_militaris-GEM\ComplementaryData\cmtKEGGModel.mat'
load 'C:\Users\Nachonase\Documents\GitHub\Cordyceps_militaris-GEM\ComplementaryData\cmtBlastYeast.mat'
cmtDraft = getModelFromHomology({tmpModel},cmtBlastYeast,'cmt',{},2,false,10^-5,100);
cmtDraft.id = 'cmtDraft';
cmtDraft.description = 'draft model of Cordyceps reconstucted by yeastHomology';
[newConnected, cannotConnect, addedRxns, cmtModel, exitFlag]=...
    fillGaps(cmtDraft,{tmpModel});

% identify intersect reactions between C. militaris and yeast based on KEGGID 
cmtKEGG = unique(cmtKEGGModel.rxns);
yeastKEGG = unique(tmpModel.rxnKEGGID);
intersectKEGG = intersect(yeastKEGG,cmtKEGG);
intersectRxnsInKEGG = tmpModel.rxns(find(ismember(tmpModel.rxnKEGGID,intersectKEGG)));
add_intersectRxnsInKEGG = setdiff(intersectRxnsInKEGG,cmtModel.rxns);
cmtModel = addRxnsGenesMets(cmtModel,tmpModel,add_intersectRxnsInKEGG,false,...
    'intersect reactions between C. militaris and yeast based on KEGGID',4); 

% merge core into draft 
cmtModel = addRxnsGenesMets(cmtModel,cmtCore,cmtCore.rxns,false,...
    'essential reactions for C. militaris growth',1); 

cmtModel.equations = constructEquations(cmtModel);

id_cmt_ex_c = regexp(cmtModel.rxns,...
    'cmt_ex_c_','all');
cmt_ex_c = cmtModel.rxns(find(~cellfun('isempty',id_cmt_ex_c)));
cmtModel=setParam(cmtModel,'lb',cmt_ex_c,0);

cmtModel=setParam(cmtModel,'lb','cmt_ex_c_glucose',-1);
cmtModel=setParam(cmtModel,'obj','cmt_ex_growth',1);
sol = solveLP(cmtModel);
cmtModel.haveFlux(1) = sol;
for i = 1:numel(cmt_ex_c)
    cmtModel=setParam(cmtModel,'lb',cmt_ex_c,0);
    cmtModel=setParam(cmtModel,'lb',cmt_ex_c(i),-1);
    sol = solveLP(cmtModel);
    cmtModel.haveFlux(i) = sol;
end

cmtModel.compareFlux = table(cmtModel.rxns,cmtModel.equations,cmtModel.haveFlux(1).x,...
    cmtModel.haveFlux(2).x,cmtModel.haveFlux(3).x,cmtModel.haveFlux(4).x,...
    cmtModel.haveFlux(5).x,cmtModel.haveFlux(6).x,cmtModel.haveFlux(7).x,...
    cmtModel.haveFlux(8).x,cmtModel.haveFlux(9).x,cmtModel.haveFlux(10).x,...
    cmtModel.haveFlux(11).x);

cmtModel=setParam(cmtModel,'lb',cmt_ex_c,0);
cmtModel=setParam(cmtModel,'lb','cmt_ex_c_glucose',-1);
sol = solveLP(cmtModel);
printFluxes(cmtModel,sol.x,true)

cmtModel.id = 'cmtModel';
cmtModel.description = 'functional cmtModel with cmtBiomass';

%% stage 4: call_mileStone44 > assign comps.
%compartment extention
growthRxn = {'r_2111'};
exchangeRxns = setdiff(getExchangeRxns(yeastGEM),growthRxn);
transportRxnsIdx = getTransportRxns(yeastGEM);
transportRxns = setdiff(yeastGEM.rxns(transportRxnsIdx),growthRxn);
non_comps = union(exchangeRxns,transportRxns);
non_comps = union(non_comps,growthRxn); % all unassign compartment rxns in Class I
% There are 4 groups depend on their compartments.
cytoplasmRxns = setdiff(yeastGEM.rxns(getRxnsInComp(yeastGEM,'c',false)),non_comps);
extracellularRxns = setdiff(yeastGEM.rxns(getRxnsInComp(yeastGEM,'e',false)),non_comps);
mitochondrionRxns = setdiff(yeastGEM.rxns(getRxnsInComp(yeastGEM,'m',false)),non_comps);
peroxisomeRxns = setdiff(yeastGEM.rxns(getRxnsInComp(yeastGEM,'p',false)),non_comps);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tracking their group names into their rxn IDs
model = yeastGEM; %track compartment infomation into rxn ID
% considering rxns in Class I, first!
for i = 1:numel(model.rxns)
    if ismember(model.rxns{i},growthRxn);
        model.rxns{i} = strcat(model.rxns{i},'-growthRxn'); % this growthRxn is an exchage rxn which occurs in cytoplasm.
    else
        if ismember(model.rxns{i},exchangeRxns);
        model.rxns{i} = strcat(model.rxns{i},'-exchangeRxns'); %if rxn is an exchage rxn but not growth are in extracellular
        else
            if ismember(model.rxns{i},transportRxns);
                model.rxns{i} = strcat(model.rxns{i},'-transportRxns'); %if rxn can not assign any compartment
            else
                % if the rxn is not found in any groups of Class I, 
                % it is a general rxn which located in certain compartments!
                if ismember(model.rxns{i},mitochondrionRxns);
                    model.rxns{i} = strcat(model.rxns{i},'-mitochondrionRxns'); %if rxn located in mitochondrion
                else
                    if ismember(model.rxns{i},peroxisomeRxns);
                        model.rxns{i} = strcat(model.rxns{i},'-peroxisomeRxns'); %if rxn located in peroxisome
                    else
                        if ismember(model.rxns{i},extracellularRxns);
                            model.rxns{i} = strcat(model.rxns{i},'-extracellularRxns'); 
                        else
                            % The rxns that are not included in any groups 
                            % will be assigned into cytoplasm.
                            model.rxns{i} = strcat(model.rxns{i},'-cytoplasmRxns'); 
                        end
                    end
                end
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% The tracking group names in the new rxn IDs will use after mergeCompartments
[model1, deletedRxns, duplicateRxns] = mergeCompartments(model); % must keep duplicateRxns
model2 = mergeCompartments(yeastGEM); 
model1.allRxns = cell(numel(model1.rxns),1);
model1.allRxns = strcat(model1.rxns,';',duplicateRxns); % Some rxns found in mergeCompartments model
                                                        % have duplicate rxns in different compartments. 
                                                        % Concatenated RxnID with their duplicateRxns
                                                        % has been used for
                                                        % re-grouping each
                                                        % rxn into each
                                                        % compartment
                                                        % again.
                                                        

rxnsTo_c_Idx = regexp(model1.allRxns,...
    'cytoplasmRxns','all');
rxnsTo_c_Idx = find(~cellfun('isempty',rxnsTo_c_Idx));
rxnsTo_c = model2.rxns(rxnsTo_c_Idx); % re-group reactions in cytoplasm

rxnsTo_m_Idx = regexp(model1.allRxns,...
    'mitochondrionRxns','all');
rxnsTo_m_Idx = find(~cellfun('isempty',rxnsTo_m_Idx));
rxnsTo_m = model2.rxns(rxnsTo_m_Idx); % re-group reactions in mitochondria

rxnsTo_p_Idx = regexp(model1.allRxns,...
    'peroxisomeRxns','all');
rxnsTo_p_Idx = find(~cellfun('isempty',rxnsTo_p_Idx)); 
rxnsTo_p = model2.rxns(rxnsTo_p_Idx); % re-group reactions in peroxisome

rxnsTo_e_Idx = regexp(model1.allRxns,...
    'extracellularRxns','all');
rxnsTo_e_Idx = find(~cellfun('isempty',rxnsTo_e_Idx)); 
rxnsTo_e = model2.rxns(rxnsTo_e_Idx); % re-group reactions in extracellular

rxnsTo_t_Idx = regexp(model1.allRxns,...
    'transportRxns','all');
rxnsTo_t_Idx = find(~cellfun('isempty',rxnsTo_t_Idx));
rxnsTo_t = model2.rxns(rxnsTo_t_Idx); % re-group reactions in transportRxns

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Grouping rxns in mileStone model into 5 groups:
mileStone4 = cmtModel;
% copyTo_2c: a group of reactions in mileStone that intersect with
%            reactions in cytoplasmRxns group and growthRxn.
copyTo_c1 = union(intersect(cmtModel.rxns,rxnsTo_c),intersect(cmtModel.rxns,growthRxn));
copyTo_c2 = setdiff(cmtModel.rxns,model2.rxns); %cordycepin and mannitol reactions are also included in this group.
copyTo_2c = union(copyTo_c1,copyTo_c2);


copyTo_e = union(intersect(cmtModel.rxns,rxnsTo_e),intersect(cmtModel.rxns,exchangeRxns));
copyTo_m = intersect(cmtModel.rxns,rxnsTo_m);
copyTo_p = intersect(cmtModel.rxns,rxnsTo_p);

% copyTo_t: This group of reactions in mileStone are transportRxns which
%           must replace by their original rxns from yeastGEM_v8.0.2.
copyTo_t = intersect(cmtModel.rxns,rxnsTo_t);

mileStone4 = copyToComps(mileStone4,'c',copyTo_2c,false,{'cytoplasm'});
mileStone4 = copyToComps(mileStone4,'e',copyTo_e,false,{'extracellular'});
mileStone4 = copyToComps(mileStone4,'m',copyTo_m,false,{'mitochondrion'});                
mileStone4 = copyToComps(mileStone4,'p',copyTo_p,false,{'peroxisome'});
mileStone4 = copyToComps(mileStone4,'t',copyTo_t,false,{'transport'}); % WARNING!

mileStone44 = removeReactions(mileStone4,cmtModel.rxns,true,true,true); % remove all reactions in 's' compartment!


id_cmt_ex = regexp(mileStone44.rxns,...
    'cmt_ex_','all');
cmt_ex = mileStone44.rxns(find(~cellfun('isempty',id_cmt_ex)));
moveRxns = setdiff(cmt_ex,'cmt_ex_growth_c');
mileStone44=copyToComps(mileStone44,'e',moveRxns,false);
mileStone44 = removeReactions(mileStone44,moveRxns,true,true,true);
% Remove old transport rxns in mileStone44 (rxnsTo_t)\because they are duplicated and incorrect.
mileStone44 = removeReactions(mileStone44,mileStone44.rxns(getRxnsInComp(mileStone44,'t',true)),...
    true,true,true);
mustHave = {'r_1126';'r_1138';'r_1127';'r_1265';'r_1237';'r_1099';'r_0507';'r_0508';'r_0770';'r_0439';'r_1110';'r_0226';'r_0438';'r_1118';'r_0713';'r_0530';'r_1112';'r_1128';'r_1081';'r_4046';'r_1083';'r_1861';'r_1652'};
mileStone44 = addRxnsGenesMets(mileStone44,yeastGEM,mustHave,false,...
    'NGAM rxns yeast model',1);
mileStone44.equations = constructEquations(mileStone44);

%% stage 5.1: call_cmtCompartments > automatic add transport

[model, paddedRxns]=addTransport(mileStone44,'p','c');
[model, maddedRxns]=addTransport(model,'m','c');
[model, eaddedRxns]=addTransport(model,'e','c');
model.equations = constructEquations(model);
noneed_transport = {'T_p_to_c_s_0045_p';'T_p_to_c_s_0051_p';'T_p_to_c_s_0243_p';'T_p_to_c_s_0250_p';'T_p_to_c_s_0253_p';'T_p_to_c_s_0257_p';'T_p_to_c_s_0282_p';'T_p_to_c_s_0677_p';'T_p_to_c_s_0779_p';'T_p_to_c_s_1616_p';'T_p_to_c_s_1620_p';'T_p_to_c_s_2792_p';'T_p_to_c_s_2793_p';'T_p_to_c_s_2794_p';'T_p_to_c_s_2795_p';'T_p_to_c_s_2796_p';'T_p_to_c_s_2812_p';'T_p_to_c_s_2813_p';'T_m_to_c_s_0033_m';'T_m_to_c_s_0116_m';'T_m_to_c_s_0162_m';'T_m_to_c_s_0166_m';'T_m_to_c_s_0169_m';'T_m_to_c_s_0185_m';'T_m_to_c_s_0190_m';'T_m_to_c_s_0197_m';'T_m_to_c_s_0215_m';'T_m_to_c_s_0234_m';'T_m_to_c_s_0274_m';'T_m_to_c_s_0282_m';'T_m_to_c_s_0286_m';'T_m_to_c_s_0315_m';'T_m_to_c_s_0428_m';'T_m_to_c_s_0430_m';'T_m_to_c_s_0432_m';'T_m_to_c_s_0677_m';'T_m_to_c_s_0680_m';'T_m_to_c_s_0709_m';'T_m_to_c_s_0710_m';'T_m_to_c_s_0714_m';'T_m_to_c_s_0734_m';'T_m_to_c_s_0748_m';'T_m_to_c_s_0830_m';'T_m_to_c_s_0832_m';'T_m_to_c_s_0847_m';'T_m_to_c_s_0850_m';'T_m_to_c_s_0853_m';'T_m_to_c_s_0924_m';'T_m_to_c_s_0929_m';'T_m_to_c_s_0932_m';'T_m_to_c_s_0937_m';'T_m_to_c_s_0943_m';'T_m_to_c_s_0965_m';'T_m_to_c_s_0969_m';'T_m_to_c_s_1006_m';'T_m_to_c_s_1025_m';'T_m_to_c_s_1029_m';'T_m_to_c_s_1032_m';'T_m_to_c_s_1048_m';'T_m_to_c_s_1051_m';'T_m_to_c_s_1077_m';'T_m_to_c_s_1099_m';'T_m_to_c_s_1101_m';'T_m_to_c_s_1148_m';'T_m_to_c_s_1219_m';'T_m_to_c_s_1222_m';'T_m_to_c_s_1311_m';'T_m_to_c_s_1314_m';'T_m_to_c_s_1318_m';'T_m_to_c_s_1384_m';'T_m_to_c_s_1386_m';'T_m_to_c_s_1403_m';'T_m_to_c_s_1405_m';'T_m_to_c_s_1413_m';'T_m_to_c_s_1416_m';'T_m_to_c_s_1491_m';'T_m_to_c_s_1527_m';'T_m_to_c_s_1529_m';'T_m_to_c_s_1533_m';'T_m_to_c_s_1535_m';'T_m_to_c_s_1537_m';'T_m_to_c_s_1561_m';'T_m_to_c_s_1583_m';'T_m_to_c_s_1585_m';'T_m_to_c_s_1587_m';'T_m_to_c_s_1591_m';'T_m_to_c_s_1594_m';'T_m_to_c_s_1596_m';'T_m_to_c_s_1598_m';'T_m_to_c_s_1600_m';'T_m_to_c_s_1602_m';'T_m_to_c_s_1604_m';'T_m_to_c_s_1608_m';'T_m_to_c_s_1610_m';'T_m_to_c_s_1612_m';'T_m_to_c_s_1614_m';'T_m_to_c_s_0149_m';'T_m_to_c_s_0319_m';'T_m_to_c_s_0347_m';'T_e_to_c_s_0714_e';'T_e_to_c_s_1405_e';'T_e_to_c_s_1497_e';'T_e_to_c_s_1489_e';'T_e_to_c_s_0925'};
model = removeReactions(model,noneed_transport,true,true,true);
model.equations = constructEquations(model);

id_cmt_ex_c = regexp(model.rxns,...
    'cmt_ex_c_','all');
cmt_ex_c = model.rxns(find(~cellfun('isempty',id_cmt_ex_c)));

model=setParam(model,'lb',cmt_ex_c,0);
model=setParam(model,'lb','cmt_ex_c_glucose_c_e',-1);
model=setParam(model,'obj','cmt_ex_growth_c',1);
sol = solveLP(model);
printFluxes(model,sol.x,true);


model=setParam(model,'obj','cmt_ex_growth_c',1);
sol = solveLP(model);
cmtCompartments = model;

%% stage 5.2: call_correctedTransport > correct transport by adding for template
% Transport rxns submodel from yeastGEM_v8.0.2
transportRxnsIdx = getTransportRxns(yeastGEM);
transportRxns = yeastGEM.rxns(transportRxnsIdx);
nontransportRxns = setdiff(yeastGEM.rxns,transportRxns);
transportM = removeReactions(yeastGEM,nontransportRxns,true,true,true);
transportM.equations = constructEquations(transportM);

% get Transport rxns from all 4 compartments
transportTo_c_Idx = getRxnsInComp(transportM,'c',true);
transportTo_e_Idx = getRxnsInComp(transportM,'e',true);
transportTo_m_Idx = getRxnsInComp(transportM,'m',true);
transportTo_p_Idx = getRxnsInComp(transportM,'p',true);
transportTo_e = transportM.rxns(transportTo_e_Idx);
transportTo_c = transportM.rxns(transportTo_c_Idx);
transportTo_m = transportM.rxns(transportTo_m_Idx);
transportTo_p = transportM.rxns(transportTo_p_Idx);
%
transportBetween_ec = intersect(transportTo_e,transportTo_c);
transportBetween_mc = intersect(transportTo_m,transportTo_c);
transportBetween_pc = intersect(transportTo_p,transportTo_c);
emp = [transportBetween_ec;transportBetween_mc;transportBetween_pc];
nonEMP = setdiff(transportM.rxns,emp);
transportEMP = removeReactions(transportM,nonEMP,true,true,true);

% remove fake ATP transport
fake_nucleo = {'T_m_to_c_s_0394_m';'T_m_to_c_s_0739_m';'T_m_to_c_s_0785_m';...
    'T_m_to_c_s_0423_m';'T_m_to_c_s_0434_m';...
    'T_p_to_c_s_0423_p';'T_p_to_c_s_0434_p';...
    'T_m_to_c_s_0687_m';'T_m_to_c_s_1198_m'}

cmtGEM4 = removeReactions(cmtCompartments,fake_nucleo,true,true,true);
tran_ATPin = {'r_1111';'r_1116';'r_1276';'r_1231';'r_1131';'r_1130';...
    'r_1230';'r_1232';'r_1275';'r_1175'};
cmtGEM4 = addRxnsGenesMets(cmtGEM4,transportEMP,tran_ATPin,false,'corrected transport rxns',1);
% correct carnitine shuttle
cmtGEM4 = removeReactions(cmtGEM4,{'T_p_to_c_s_0021_p';'T_p_to_c_s_1235_p';...
    'T_m_to_c_s_0021_m';'T_m_to_c_s_1235_m'},true,true,true);
trans_carnitine = {'r_1120';'r_1191';'r_1638';'r_1673';'r_1674';'r_1882';'r_1976'};
cmtGEM4 = addRxnsGenesMets(cmtGEM4,transportEMP,trans_carnitine,false,'corrected transport rxns',1);

aaTransport = {'r_1183';'r_1184';'r_1186';'r_1190';'r_1192';'r_1196';'r_1199';'r_1173';'r_1201';'r_1211';'r_1213';'r_1214';'r_1215';'r_1216';'r_1217';'r_1218';'r_1219';'r_1223';'r_1224';'r_1205'};
cmtGEM4 = addRxnsGenesMets(cmtGEM4,transportEMP,aaTransport,false,'aa transport rxns',1);

% correct malate citrate succinate oxoloacetate shuttle
remove_malateShuttle = {'T_p_to_c_s_0066_p';'T_p_to_c_s_0180_p';...
    'T_p_to_c_s_0522_p';'T_p_to_c_s_0940_p';'T_p_to_c_s_1271_p';...
    'T_m_to_c_s_0066_m';'T_m_to_c_s_0180_m';'T_m_to_c_s_0522_m';...
    'T_m_to_c_s_0940_m';'T_m_to_c_s_1271_m';'T_m_to_c_s_1458_m';...
    'T_e_to_c_s_0181';'T_p_to_c_s_0633_p';'T_m_to_c_s_0516_m';...
    'T_m_to_c_s_0633_m';'T_m_to_c_s_0725_m';'T_m_to_c_s_0834_m';...
    'T_m_to_c_s_1322_m';'T_m_to_c_s_1360_m';'T_m_to_c_s_1399_m'};
tran_malate = {'r_2132';'r_1689';'r_1226';'r_1930';'r_1688';...
    'r_1239';'r_1245';'r_1264'};
cmtGEM4 = removeReactions(cmtGEM4,remove_malateShuttle,true,true,true);
cmtGEM4 = addRxnsGenesMets(cmtGEM4,transportEMP,tran_malate,false,'corrected transport rxns',1);
remove_aspartateShuttle = {'T_p_to_c_s_0973_p';'T_p_to_c_s_0991_p'};
cmtGEM4 = removeReactions(cmtGEM4,remove_aspartateShuttle,true,true,true);
cmtGEM4 = addRxnsGenesMets(cmtGEM4,transportEMP,{'r_1659';'r_1657';'r_2034'},false,'corrected transport rxns',1);

remove_pCoA = {'T_p_to_c_s_0054_p';'T_p_to_c_s_0229_p';'T_p_to_c_s_0367_p';'T_p_to_c_s_0373_p';
'T_p_to_c_s_1073_p';'T_p_to_c_s_1176_p';'T_p_to_c_s_1262_p';'T_p_to_c_s_1302_p';
'T_p_to_c_s_1454_p';'T_p_to_c_s_1479_p';'T_p_to_c_s_1513_p';'T_p_to_c_s_1516_p';
'T_p_to_c_s_1519_p';'T_p_to_c_s_2814_p';'T_p_to_c_s_2819_p';'T_p_to_c_s_0816_p'};
cmtGEM4 = removeReactions(cmtGEM4,remove_pCoA);

add_newTransport = {'r_1134';'r_1135';'r_1139';'r_1166';'r_1227';'r_1244';'r_2079';'r_1258';'r_1707';'r_1717';'r_1719';'r_1115';'r_1266';'r_1277';'r_1697';'r_1824';'r_1979';};
remove_oldTransport = {'T_e_to_c_s_0419_e';'T_e_to_c_s_0445_e';'T_e_to_c_s_0456_e';'T_e_to_c_s_0553_e';'T_e_to_c_s_0563_e';'T_e_to_c_s_0793_e';'T_e_to_c_s_0803_e';'T_e_to_c_s_1275_e';'T_e_to_c_s_1322_e';'T_e_to_c_s_1467_e';'T_e_to_c_s_0561_c_e';'T_e_to_c_s_0571_c_e';'T_e_to_c_s_1105_c_e';'T_e_to_c_s_1520_c_e';'T_e_to_c_s_0578_c_e';'T_e_to_c_s_0548_c_e';'T_e_to_c_s_0067'};
cmtGEM4 = removeReactions(cmtGEM4,remove_oldTransport);
cmtGEM4 = addRxnsGenesMets(cmtGEM4,transportEMP,add_newTransport,false,'corrected transport rxns',1);

indexFructose = find(ismember(cmtGEM4.rxns,'r_1134'));
cmtGEM4.S(find(ismember(cmtGEM4.mets,'s_0793_e')),indexFructose) = 0;
cmtGEM4.S(find(ismember(cmtGEM4.mets,'s_0793_c')),indexFructose) = 0;
        
noused = {'T_m_to_c_s_0750_m';'T_p_to_c_s_0204_p';'T_p_to_c_s_1051_p';...
    'T_m_to_c_s_0025_m';'T_m_to_c_s_0973_m';'T_m_to_c_s_0991_m';...
    'T_m_to_c_s_0176_m';'T_m_to_c_s_0304_m';'T_m_to_c_s_0793_m';...
    'T_m_to_c_s_1212_m';'T_m_to_c_s_1487_m';'T_m_to_c_s_0120_m';...
    'T_m_to_c_s_0306_m';'T_m_to_c_s_0803_m';'T_m_to_c_s_1003_m';...
    'T_m_to_c_s_1207_m';'T_m_to_c_s_0837_m';...
    'T_m_to_c_s_0232_m';'T_m_to_c_s_0056_m';...
    'T_m_to_c_s_0362_m';'r_1024_e';'T_e_to_c_Cellobiose_c_e';'T_e_to_c_Lactose_c_e';'T_e_to_c_Cordycepin_c_e';...
    'T_m_to_c_s_0359_m';'T_m_to_c_s_0625_m';'T_m_to_c_s_0689_m';'T_m_to_c_s_0754_m';...
    'T_m_to_c_s_1203_m';'T_p_to_c_s_0529_p';'T_p_to_c_s_1260_p';'T_p_to_c_s_1293_p'};
cmtGEM4 = removeReactions(cmtGEM4,noused);
cmtGEM4.equations = constructEquations(cmtGEM4);

in.rxns = {'T_p_to_c_s_0456_p';'T_p_to_c_s_0793_p';'T_p_to_c_s_0803_p';'T_p_to_c_s_0837_p';'T_p_to_c_s_1065_p';'T_p_to_c_s_1161_p';'T_p_to_c_s_1198_p';'T_p_to_c_s_1203_p';'T_p_to_c_s_1207_p';'T_p_to_c_s_1212_p';'T_p_to_c_s_1275_p';'T_p_to_c_s_1286_p';'T_p_to_c_s_1449_p';'T_m_to_c_s_0010_m';'T_m_to_c_s_0118_m';'T_m_to_c_s_0178_m';'T_m_to_c_s_0218_m';'T_m_to_c_s_0291_m';'T_m_to_c_s_0349_m';'T_m_to_c_s_0367_m';'T_m_to_c_s_0373_m';'T_m_to_c_s_0419_m';'T_m_to_c_s_0445_m';'T_m_to_c_s_0456_m';'T_m_to_c_s_0529_m';'T_m_to_c_s_0551_m';'T_m_to_c_s_0629_m';'T_m_to_c_s_0722_m';'T_m_to_c_s_0767_m';'T_m_to_c_s_0955_m';'T_m_to_c_s_1016_m';'T_m_to_c_s_1021_m';'T_m_to_c_s_1039_m';'T_m_to_c_s_1045_m';'T_m_to_c_s_1056_m';'T_m_to_c_s_1266_m';'T_m_to_c_s_1275_m';'T_m_to_c_s_1616_m';'T_m_to_c_s_1620_m'};
cmtGEM4=setParam(cmtGEM4,'lb',in.rxns,-1000);
cmtGEM4=setParam(cmtGEM4,'ub',in.rxns,1000);
more.rxns = {'cordycepinTransport';'sucroseTransport';'sucrase';...
    'CellobioseTransport';'LactoseTransport';'gabaPermease';'nicotinatePermease'};
more.equations = {'Cordycepin[c] => Cordycepin[e]';...
    'H+[e] + sucrose[e] => H+[c] + sucrose[c]';...
    'H2O[c] + sucrose[c] => D-fructose[c] + D-glucose[c]';...
    'H+[e] + Cellobiose[e] => H+[c] + Cellobiose[c]';...
    'H+[e] + Lactose[e] => H+[c] + Lactose[c]';...
    'H+[e] + gamma-aminobutyrate[e] => H+[c] + gamma-aminobutyrate[c]';...
    'H+[e] + nicotinate[e] => H+[c] + nicotinate[c]'};
cmtGEM4 = addRxns(cmtGEM4,more,3,'',true);
cmtGEM4=setParam(cmtGEM4,'lb',more.rxns,0);
cmtGEM4=setParam(cmtGEM4,'ub',more.rxns,1000);
cmtGEM4.equations = constructEquations(cmtGEM4);


%% stage 6: Fit GAM and NGAM
id_cmt_ex_c = regexp(cmtGEM4.rxns,'cmt_ex_c_','all');
cmt_ex_c = cmtGEM4.rxns(find(~cellfun('isempty',id_cmt_ex_c)));
cmtGEM4=setParam(cmtGEM4,'lb',cmt_ex_c,0);
cmtGEM4=setParam(cmtGEM4,'lb','r_4046',0);
cmtGEM4=setParam(cmtGEM4,'ub','r_4046',1000);
cmtGEM4=setParam(cmtGEM4,'obj','r_4046',1);
cmtGEM4=setParam(cmtGEM4,'lb','cmt_ex_c_glucose_c_e',-0.02);

cmtGEM4=setParam(cmtGEM4,'eq','cmt_ex_growth_c',0);
sol = solveLP(cmtGEM4);
cmtNGAM=setParam(cmtGEM4,'eq','r_4046',-sol.f);
cmtNGAM=setParam(cmtNGAM,'lb','r_1654_e',-1.4);
cmtNGAM=setParam(cmtNGAM,'lb','cmt_ex_c_glucose_c_e',-0.15);
cmtNGAM=setParam(cmtNGAM,'lb','cmt_ex_growth_c',0);
cmtNGAM=setParam(cmtNGAM,'ub','cmt_ex_growth_c',1000);
cmtNGAM=setParam(cmtNGAM,'obj','cmt_ex_growth_c',1);
solveLP(cmtNGAM);

%% stage 7: add genes
[~,grRules] = xlsread('C:\Users\Nachonase\Documents\GitHub\Cordyceps_militaris-GEM\ComplementaryData\genesToAdd.xlsx','grRules'); 
[~,genes] = xlsread('C:\Users\Nachonase\Documents\GitHub\Cordyceps_militaris-GEM\ComplementaryData\genesToAdd.xlsx','genesToAdd'); 
genesToAdd = struct();
genesToAdd.genes = genes(1:end,1);
grRulesToAdd.rxns = grRules(1:end,1);
grRulesToAdd.grRules = grRules(1:end,2);

cmtFinal = addGenesRaven(cmtNGAM,genesToAdd);
[a, b] = ismember(cmtFinal.rxns,grRulesToAdd.rxns);
I = find(a);
cmtFinal.grRules(I) = grRulesToAdd.grRules(b(I));
cmtFinal.mets{1218} = 'GABA_e';
cmtFinal.mets{1219} = 'nicotinate_e';
cmtFinal.mets{1220} = 'sucrose_c';


cmtForGit = cmtFinal;
cmtForGit.subSystems = cell(numel(cmtForGit.subSystems),1);
cmtForGit.rxnMiriams = cell(numel(cmtForGit.rxnMiriams),1);
exportForGit(cmtForGit,'model','C:\Users\Nachonase\Documents\GitHub\Cordyceps_militaris-GEM\');

%% end %%
