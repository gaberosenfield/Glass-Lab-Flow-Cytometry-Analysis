% This is the script Pedro Goncalves used to analyze flow cytometry data
% in his 2019 paper.  It requires the export_fig functions from the matlab
% file exchange if you want to output example figures.
function FC_analysis_Goncalves(makehistplots)
% This script removes ungerminated conidia from the matched germinated samples in FCS data.
% It then uses Otsu's method to automatically gate the fluorescence data for each sample independently.
% Finally, it corrects the percentage of highly fluorescent cells based on exponential decay fitting of data
% above the fluorescence gate.
close all;

% Set Control Variables below
conPrctile = 95; % this variable controls the forward scatter area threshold that differentiates conidia from germlings
makeFigs = false; % change this variable to true to output gating example figures
figLoc = '/Users/Gabe/Dropbox/UC Berkeley/Glass Lab/Literature/My Stuff/Goncalves et al 2019/Figures/'; %change this variable to the directory path into which you want example figures to be saved

% Validate input
if nargin
    if makehistplots % this argument controls whether plots will be generated.
        makehistplots = true;
    else
        makehistplots = false;
    end
else
    makehistplots = false; % do not make histogram plots by default.
end

fcsfiles = dir('*.fcs');
Nfcsfiles = length(fcsfiles);
crossfiles_sts = cell(Nfcsfiles,1);
conidiafiles_sts = cell(Nfcsfiles,2);
strainfiles_sts = cell(Nfcsfiles,1);
icross=0;
iconidia=0;
istrain = 0;

for ifile = 1:Nfcsfiles
    fcsfile_ifile = fcsfiles(ifile).name;
    match_sts=regexp(fcsfile_ifile,' C_\d\d\d.fcs', 'once');
    if isempty(match_sts)
        match_sts = regexp(fcsfile_ifile,' x ', 'once');
        if isempty(match_sts)
            istrain = istrain+1;
            strainfiles_sts{istrain} = fcsfiles(ifile).name;
        else
            icross = icross+1;
            crossfiles_sts{icross} = fcsfiles(ifile).name;
        end
    else
        iconidia = iconidia+1;
        conidiafiles_sts{iconidia,1} = fcsfiles(ifile).name;
    end
end

conidiafiles_sts(iconidia+1:end,:)=[];
crossfiles_sts(icross+1:end)=[];
strainfiles_sts(istrain+1:end) = [];

N_conidia = iconidia;
N_crosses = icross;
N_strain = istrain;

% Set thresholds on Forward Scatter for something to be counted as a germling (i.e. not as a conidia)
for iconidia = 1:N_conidia
    [conidiadat,~] = fca_readfcs(conidiafiles_sts{iconidia,1});
    conidia_FSC_sort = sort(conidiadat(:,1),'descend');
    N_events = length(conidiadat(:,1));
    FSCthres  = conidia_FSC_sort(round(N_events*(1-conPrctile/100))); % determine the forward scatter area threshold that differentiates conidia from germlings
    conidiafiles_sts{iconidia,2} = FSCthres;
    
    if makeFigs % separately record the conidial FSCA data, SSCA data, & gate from 2489 & S9S
        if contains(conidiafiles_sts{iconidia,1},'2489 C_')
            c2489FSCA = conidiadat(:,1);
            c2489SSCA = conidiadat(:,3);
            FSCthresWT = FSCthres;
            conidia_SSC_sort = sort(conidiadat(:,3),'descend');
            SSCthres2489 = conidia_SSC_sort(round(N_events*(1-conPrctile/100)));
        elseif contains(conidiafiles_sts{iconidia,1},'S9S C_')
            cS9SFSCA = conidiadat(:,1);
            FSCthresS9S = FSCthres;
        end
    end
end

outputFID = fopen('strain specific threshold w scoring.csv','w');

fprintf(outputFID,['Strain 1,Strain 2,# germinated,SB threshold,PI threshold,# dead by SB,# dead by PI,' ...
    '# dead by both,%% dead by SB,%% dead by PI,%% dead by both,SB error score,PI error score,' ...
    'corrected %% dead by SB,corrected %% dead by PI,\n']);

% Loop through the individual strain samples and analyze the data
for istrain = 1:N_strain
    fcs_st = strainfiles_sts{istrain};
    match_undersc = regexp(fcs_st,'_');
    strain_st = fcs_st(match_undersc(1)+1:match_undersc(end)-1);
    
    FSCthres = [];
    for i_conidia=1:N_conidia
        ismatch = regexp(conidiafiles_sts{i_conidia,1},strcat("_",strain_st), 'once');
        if ~isempty(ismatch)
            FSCthres = conidiafiles_sts{i_conidia,2};
        end
    end
    
    if isempty(FSCthres)
       error(['I could not find the file containing conidial data for ' strain_st]);
    end

    % The if statement below generates figures demonstrating the conidial gating method used in this program
    if makeFigs && strcmp(strain_st,'2489')
        [wtDat,~] = fca_readfcs(strainfiles_sts{istrain});
        wtFSCA = wtDat(:,1);
        germFSCA = wtFSCA(wtFSCA>FSCthresWT);
        wtSSCA = wtDat(:,3);
        aboveFSCgate = 100*length(germFSCA)/length(wtFSCA);
        aboveSSCgate = 100*sum(wtSSCA>SSCthres2489)/length(wtSSCA);
        aboveBothGates = 100*length(intersect(find(wtFSCA>FSCthresWT),find(wtSSCA>SSCthres2489)))/length(wtFSCA);
        
        figA = figure('Units','inches','Position',[0 0 8.5 8.5]); % Generate figure showing why we don't use SSCA for gating conidia
        scatter(log(wtFSCA),log(wtSSCA),'.','MarkerEdgeColor',[0 0.4470 0.7410]);
        hold on;
        scatter(log(c2489FSCA),log(c2489SSCA),'.','MarkerEdgeColor',[0.8500 0.3250 0.0980]);
        axis manual;
        xl = xlim;
        yl = ylim;
        lFit = polyfit(log(wtFSCA),log(wtSSCA),1);
        regLine = polyval(lFit,log(wtFSCA));
        plot(log(wtFSCA),regLine,'-','LineWidth',2,'Color',[0.4660 0.6740 0.1880]);
        cc = corrcoef(log(wtFSCA),log(wtSSCA));
        plot(log(FSCthresWT)*[1 1],[yl(1) yl(2)],'--k','LineWidth',2);
        plot([xl(1) xl(2)],log(SSCthres2489)*[1 1],'--','LineWidth',2,'Color',[0.4940 0.1840 0.5560]);
        xlabel('Natural Log of Forward Scatter Area (arbitrary units)','FontSize',15);
        ylabel('Natural Log of Side Scatter Area (arbitrary units)','FontSize',15);
        lg = legend({['FGSC2489 cells, n = ' num2str(length(wtFSCA))],['FGSC2489 conidia, n = ' num2str(length(c2489FSCA))], ...
            ['Linear regression of FGSC2489 4h tpi FSC-A vs SSC-A, R^{2} = ' num2str(cc(1,2),2)], ...
            [num2str(aboveFSCgate,2) '% of FGSC2489 4h tpi above ' num2str(conPrctile) 'th percentile conidial FSC-A gate'], ...
            [num2str(aboveSSCgate,2) '% of FGSC2489 4h tpi above ' num2str(conPrctile) 'th percentile conidial SSC-A gate']});
        lg.Location = 'northwest';
        lg.FontSize = 12;
        title({'FGSC2489 FSC-A vs SSC-A Scatter Plot\rm',[num2str(aboveBothGates,2) ...
            '% of FGSC2489 4h tpi above both gates']},'FontSize',15);
        figA.Color = [1 1 1];
        export_fig(figA,[figLoc 'FGSC2489 FSC-A vs SSC-A'],'-png','-jpg');
        print(figA,'-dpdf',[figLoc 'FGSC2489 FSC-A vs SSC-A']);
        print(figA,'-dsvg',[figLoc 'FGSC2489 FSC-A vs SSC-A']);
        close all;
        
        figB = figure('Units','inches','Position',[0 0 8.5 8.5]); % Generate figure showing how conidia are gated from a single strain
        wth = histogram(log(wtFSCA),200,'EdgeAlpha',0);
        hold on;
        histogram(log(c2489FSCA),'BinEdges',wth.BinEdges,'LineWidth',2,'EdgeAlpha',0,'FaceAlpha',0.6);
        axis manual;
        yl = ylim;
        plot(log(FSCthresWT)*[1 1],[yl(1) yl(2)],'--k','LineWidth',2);
        histogram(log(germFSCA),'BinEdges',wth.BinEdges,'DisplayStyle','stairs','LineWidth',1,'EdgeColor',[0.4940 0.1840 0.5560]);
        xlabel('Natural Log of Forward Scatter Area (arbitrary units)','FontSize',15);
        ylabel('Counts','FontSize',15);
        lg = legend({['FGSC2489 4h tpi FSC-A, n = ' num2str(length(wtFSCA))],['FGSC2489 conidia FSC-A, n = ' ...
            num2str(length(c2489FSCA))],[num2str(conPrctile) 'th percentile of conidial FSC-A distribution'], ...
            [num2str(aboveFSCgate,2) '% of FGSC2489 4h tpi above ' newline num2str(conPrctile) 'th percentile conidial FSC-A gate']});
        lg.FontSize = 12;
        lg.Position = [0.4314 0.7500 0.4608 0.1600];
        drawnow;
        lg.EntryContainer.NodeChildren(1).Icon.Transform.Children.Children(1).VertexData(1:2,:) = [0 1 0 1;0.5 0.5 0.5 0.5];
        title({'Conidial Gating Example: FGSC2489 FSC-A Histograms\rm', ...
            [num2str(aboveFSCgate,2) '% of FGSC2489 4h tpi germinated']},'FontSize',15);
        figB.Color = [1 1 1];
        export_fig(figB,[figLoc 'FGSC2489 Conidial Gating'],'-png','-jpg');
        print(figB,'-dpdf',[figLoc 'FGSC2489 Conidial Gating']);
        print(figB,'-dsvg',[figLoc 'FGSC2489 Conidial Gating']);
        close all;
    end
    
    findNdead(fcs_st,FSCthres,outputFID,strain_st,'',makehistplots,makeFigs,figLoc); % analyze the data from this iteration
end

% Loop through the mixed samples and analyze the data
for icross = 1:N_crosses
    fcs_st = crossfiles_sts{icross};
    match_undersc = regexp(fcs_st,'_'); % starting underscore
    match_x = regexp(fcs_st,' x '); % starting underscore
    strain1_st = fcs_st(match_undersc(1)+1:match_x-1);
    strain2_st = fcs_st(match_x+3:match_undersc(end)-1);
    
    FCSthres = 0;
    N_FCSthres = 0;
    for i_conidia=1:N_conidia
        ismatch = [regexp(conidiafiles_sts{i_conidia,1},strcat("_",strain1_st),'once') ...
            regexp(conidiafiles_sts{i_conidia,1},strcat("_",strain2_st),'once')];
        if ~isempty(ismatch)
            FSCthres = max(FCSthres,conidiafiles_sts{i_conidia,2});
            N_FCSthres = N_FCSthres+1;
        end
    end
    
    if ~any(N_FCSthres)
        error(['I could not find either of the files containing conidial data for ' strain1_st ' & ' strain2_st '.']);
    elseif N_FCSthres == 1
        error(['I could not find one of the files containing conidial data for ' strain1_st ' or ' strain2_st '.']);
    end
    
    % The if statement below generates figures demonstrating the conidial gating method used in this program
    if makeFigs && strcmp(strain1_st,'2489') && strcmp(strain2_st,'S9S')
        [mixDat,~] = fca_readfcs(crossfiles_sts{icross});
        mixFSCA = mixDat(:,1);
        germMixFSCA = mixFSCA(mixFSCA>FSCthres);
        aboveFSCgate = 100*length(germMixFSCA)/length(mixFSCA);
        
        figC = figure('Units','inches','Position',[0 0 8.5 8.5]); % Generate figure showing how conidia are gated from a mixture of strains
        mh = histogram(log(mixFSCA),200,'EdgeAlpha',0);
        hold on;
        histogram(log(c2489FSCA),'BinEdges',mh.BinEdges,'LineWidth',2,'EdgeAlpha',0,'FaceAlpha',0.6);
        histogram(log(cS9SFSCA),'BinEdges',mh.BinEdges,'LineWidth',2,'EdgeAlpha',0,'FaceAlpha',0.6);
        axis manual;
        yl = ylim;
        plot(log(FSCthresWT)*[1 1],[yl(1) yl(2)],'--k','LineWidth',2);
        plot(log(FSCthresS9S)*[1 1],[yl(1) yl(2)],'--','LineWidth',2,'Color',[0.6350 0.0780 0.1840]);
        histogram(log(germMixFSCA),'BinEdges',mh.BinEdges,'DisplayStyle','stairs','LineWidth',1,'EdgeColor',[0.4940 0.1840 0.5560]);
        xlabel('Natural Log of Forward Scatter Area (arbitrary units)','FontSize',15);
        ylabel('Counts','FontSize',15);
        lg = legend({['FGSC2489 + \itsec-9\rm swap 4h tpi FSC-A, n = ' num2str(length(mixFSCA))],['FGSC2489 conidia FSC-A, n = ' ...
            num2str(length(c2489FSCA))],['\itsec-9\rm swap conidia FSC-A, n = ' num2str(length(cS9SFSCA))], ...
            [num2str(conPrctile) 'th percentile of FGSC2489 conidial FSC-A distribution'], ...
            [num2str(conPrctile) 'th percentile of \itsec-9\rm swap conidial FSC-A distribution'], ...
            [num2str(aboveFSCgate,2) '% of FGSC2489 + \itsec-9\rm swap 4h tpi above ' newline num2str(conPrctile) 'th percentile conidial FSC-A gate']});
        lg.Location = 'northeast';
        lg.FontSize = 12;
        lg.Position = [0.3350 0.7180 0.5637 0.2000];
        drawnow;
        lg.EntryContainer.NodeChildren(1).Icon.Transform.Children.Children(1).VertexData(1:2,:) = [0 1 0 1;0.5 0.5 0.5 0.5];
        title({'Mixed Conidial Gating Example: FGSC2489 + \itsec-9\rm \bfswap FSC-A Histograms\rm', ...
            [num2str(aboveFSCgate,2) '% of FGSC2489 + \itsec-9\rm swap 4h tpi germinated']},'FontSize',15);
        figC.Color = [1 1 1];
        export_fig(figC,[figLoc 'FGSC2489 + S9S Conidial Gating'],'-png','-jpg');
        print(figC,'-dpdf',[figLoc 'FGSC2489 + S9S Conidial Gating']);
        print(figC,'-dsvg',[figLoc 'FGSC2489 + S9S Conidial Gating']);
        close all;
    end
    
    findNdead(fcs_st,FSCthres,outputFID,strain1_st,strain2_st,makehistplots,makeFigs,figLoc); % analyze the data from this iteration
end    

fclose(outputFID); %close the output file
return;

function findNdead(fcs_st,FSCthres,outputFID,strain1_st,strain2_st,makehistplots,varargin)
% Input Validation
if isempty(varargin)
    makeFigs = false; % do not make gating example figures by default
elseif length(varargin) == 1
    makeFigs = varargin{1};
    figLoc = '';
else
    makeFigs = varargin{1};
    figLoc = varargin{2};
end

[fusiondat,~] = fca_readfcs(fcs_st); % read in data
is_notconidia = fusiondat(:,1)>FSCthres; % gate out conidia

fusion_SB = fusiondat(is_notconidia,4); % fluorescence in Sytox Blue
fusion_PI = fusiondat(is_notconidia,5); % fluorescence in Propidium Iodide

N_events = sum(is_notconidia); % count # of germinated cells

is_negativesignal = fusion_SB<1 | fusion_PI<1; % remove cells with negative ln(fluorescence)
fusion_SB(is_negativesignal) = [];
fusion_PI(is_negativesignal) = [];

logfusion_SB = log(fusion_SB); % take natural log of fluorescence data
logfusion_PI = log(fusion_PI);

SBthres = multithresh_fca(logfusion_SB,1); % find SB gate (look for 1 gate separating 2 populations)
[errscore_SB,flag,sbFitVals,sbFitCounts]=otsufilter(logfusion_SB,SBthres(end)); % check to see how closely SB data above the gate approximates exponential decay
if flag==1
    warning(['Exponential fitting of SB events above threshold based on' ...
        ' bins containing <10 events for' strain1_st ' x ' strain2_st '.']);
elseif flag==2
    warning(['Not enough data for exponential fitting of SB events above' ...
        ' threshold for:' strain1_st ' x ' strain2_st '. "SB error score"' ...
        ' & "corrected %% dead by SB" will be recorded as NaN.']);
end
N_deadbySB = sum(logfusion_SB>SBthres(end)); % calculate number of germinated cells with SB fluorescence above gate

PIthres = multithresh_fca(logfusion_PI,2); % find PI gate (look for 2 gates separating 3 populations)
[errscore_PI,flag,piFitVals,piFitCounts]=otsufilter(logfusion_PI,PIthres(end)); % check to see how closely PI data above the gate approximates exponential decay
if flag==1
    warning(['Exponential fitting of PI events above threshold based on' ...
        ' bins containing <10 events for' strain1_st ' x ' strain2_st '.']);
elseif flag==2
    warning(['Not enough data for exponential fitting of PI events above' ...
        ' threshold for:' strain1_st ' x ' strain2_st '. "PI error score"' ...
        ' & "corrected %% dead by PI" will be recorded as NaN.']);
end   
N_deadbyPI = sum(logfusion_PI>PIthres(end)); % calculate number of germinated cells with PI fluorescence above gate

N_deadbyboth = sum(logfusion_SB>SBthres(end) & logfusion_PI>PIthres(end)); % calculate number of germinated cells with SB & PI fluorescence above gates

% Print output to file
fprintf(outputFID,'%s,%s,%i,%.4f,%.4f,%i,%i,%i,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,\n',strain1_st,strain2_st, ...
    N_events,SBthres(end),PIthres(end),N_deadbySB,N_deadbyPI,N_deadbyboth, ...
    (100*N_deadbySB/N_events),(100*N_deadbyPI/N_events),(100*N_deadbyboth/N_events),errscore_SB,errscore_PI, ...
    (100*(N_deadbySB/N_events)*errscore_SB),(100*(N_deadbyPI/N_events)*errscore_PI));

if makehistplots % print simple figures of each fluorescence gate to the working directory
    if isempty(strain2_st)
        plotTit = strain1_st;
    else
        plotTit = [strain1_st ' + ' strain2_st];
    end
    
    newSBthres = prctile(logfusion_SB,100-(100*(N_deadbySB/N_events)*errscore_SB));
    figure;
    histogram(logfusion_SB,200,'EdgeAlpha',0);
    hold on;
    axis manual;
    yl = ylim;
    title([plotTit ' - SB']);
    plot(SBthres(end)*[1 1],[yl(1) yl(2)],'k');
    plot(newSBthres*[1 1],[yl(1) yl(2)],'g');
    xlabel('Natural Log of Sytox Blue Fluorescence (arbitrary units)');
    ylabel('Counts');
    legend(['Germinated SB fluorescence, n = ' num2str(N_events)],'Initial SB gate','Corrected SB gate','Location','best');
    print([plotTit ' -SB.png'],'-dpng');
    close all;
    
    newPIthres = prctile(logfusion_PI,100-(100*(N_deadbyPI/N_events)*errscore_PI));
    figure;
    histogram(logfusion_PI,200,'EdgeAlpha',0,'FaceColor',[0.6350 0.0780 0.1840]);
    hold on;
    axis manual;
    yl = ylim;
    title([plotTit ' - PI']);
    plot(PIthres(end)*[1 1],[yl(1) yl(2)],'k');
    plot(newPIthres*[1 1],[yl(1) yl(2)],'g');
    xlabel('Natural Log of Propidium Iodide Fluorescence (arbitrary units)');
    ylabel('Counts');
    legend(['Germinated PI fluorescence, n = ' num2str(N_events)],'Initial PI gate','Corrected PI gate','Location','best');
    print([plotTit ' -PI.png'],'-dpng');  
    close all;   
end

if makeFigs && strcmp(strain1_st,'2489')
    if isempty(strain2_st)
        plotTit = strain1_st;
        saveTit = plotTit;
        otherTit = plotTit;
    else
        plotTit = [strain1_st ' + \itsec-9\rm \bfswap'];
        saveTit = [strain1_st ' + sec-9 swap'];
        otherTit = [strain1_st ' + \itsec-9\rm swap'];
    end
    
    % Calculate exponential decay regression of cells considered dead by SB
    logSBC = log(sbFitCounts);
    sbF = polyfit(sbFitVals,logSBC,1);
    sbPrime1 = sbF(2);
    sbPrime2 = sbF(1);
    sbDecay = exp(sbPrime1)*exp(sbFitVals*sbPrime2);
    newSBthres = prctile(logfusion_SB,100-(100*(N_deadbySB/N_events)*errscore_SB));
    
    % Calculate exponential decay regression of cells considered dead by PI
    logPIC = log(piFitCounts);
    piF = polyfit(piFitVals,logPIC,1);
    piPrime1 = piF(2);
    piPrime2 = piF(1);
    piDecay = exp(piPrime1)*exp(piFitVals*piPrime2);
    newPIthres = prctile(logfusion_PI,100-(100*(N_deadbyPI/N_events)*errscore_PI));
    
    figD = figure('Units','inches','Position',[0 0 8.5 8.5]); % Generate figure showing SB gating
    sbh = histogram(logfusion_SB,200,'EdgeAlpha',0);
    hold on;
    axis manual;
    yl = ylim;
    plot(SBthres(end)*[1 1],[yl(1) yl(2)],'--k','LineWidth',2);
    plot(sbFitVals,sbDecay,'LineWidth',2,'Color',[0.4660 0.6740 0.1880]);
    plot(newSBthres*[1 1],[yl(1) yl(2)],'--m','LineWidth',2);
    histogram(logfusion_SB(logfusion_SB>newSBthres),'BinEdges',sbh.BinEdges,'DisplayStyle','stairs','LineWidth',1,'EdgeColor',[0.4940 0.1840 0.5560]);
    xlabel('Natural Log of Sytox Blue Fluorescence (arbitrary units)','FontSize',15);
    ylabel('Counts','FontSize',15);
    lg = legend({['FGSC' otherTit ' SB, n = ' num2str(N_events)],'Initial SB gate',['Exponential decay fit of data above initial gate, R^{2} = ' ...
        num2str(1-errscore_SB,2)],'Corrected SB gate', ...
        [num2str(errscore_SB*100*N_deadbySB/N_events,2) '% of FGSC' otherTit ' above corrected SB gate']});
    lg.Location = 'northeast';
    lg.FontSize = 12;
    drawnow;
    lg.EntryContainer.NodeChildren(1).Icon.Transform.Children.Children(1).VertexData(1:2,:) = [0 1 0 1;0.5 0.5 0.5 0.5];
    title({['Sytox Blue Gating Example: FGSC' plotTit ' SB Histogram\rm'], ...
        [num2str(errscore_SB*100*N_deadbySB/N_events,2) '% of FGSC' otherTit ' dead by SB']},'FontSize',15);
    figD.Color = [1 1 1];
    export_fig(figD,[figLoc 'FGSC' saveTit ' SB Gating'],'-png','-jpg');
    print(figD,'-dpdf',[figLoc 'FGSC' saveTit ' SB Gating']);
    print(figD,'-dsvg',[figLoc 'FGSC' saveTit ' SB Gating']);
    close all;
    
    figE = figure('Units','inches','Position',[0 0 8.5 8.5]); % Generate figure showing PI gating
    pih = histogram(logfusion_PI,200,'EdgeAlpha',0,'FaceColor',[0.9000 0.0780 0.1840]);
    hold on;
    axis manual;
    yl = ylim;
    plot(PIthres(end)*[1 1],[yl(1) yl(2)],'--k','LineWidth',2);
    plot(piFitVals,piDecay,'LineWidth',2,'Color',[0.4660 0.6740 0.1880]);
    plot(newPIthres*[1 1],[yl(1) yl(2)],'--m','LineWidth',2);
    histogram(logfusion_PI(logfusion_PI>newPIthres),'BinEdges',pih.BinEdges,'DisplayStyle','stairs','LineWidth',1,'EdgeColor',[0.4940 0.1840 0.5560]);
    xlabel('Natural Log of Propidium Iodide Fluorescence (arbitrary units)','FontSize',15);
    ylabel('Counts','FontSize',15);
    lg = legend({['FGSC' otherTit ' PI, n = ' num2str(N_events)],'Initial PI gate',['Exponential decay fit of data above initial gate, R^{2} = ' ...
        num2str(1-errscore_PI,2)],'Corrected PI gate', ...
        [num2str(errscore_PI*100*N_deadbyPI/N_events,2) '% of FGSC' otherTit ' above corrected PI gate']});
    lg.Location = 'northwest';
    lg.FontSize = 12;
    drawnow;
    lg.EntryContainer.NodeChildren(1).Icon.Transform.Children.Children(1).VertexData(1:2,:) = [0 1 0 1;0.5 0.5 0.5 0.5];
    title({['Propidium Iodide Gating Example: FGSC' plotTit ' PI Histogram\rm'], ...
        [num2str(errscore_PI*100*N_deadbyPI/N_events,2) '% of FGSC' otherTit ' dead by PI']},'FontSize',15);
    figE.Color = [1 1 1];
    export_fig(figE,[figLoc 'FGSC' saveTit ' PI Gating'],'-png','-jpg');
    print(figE,'-dpdf',[figLoc 'FGSC' saveTit ' PI Gating']);
    print(figE,'-dsvg',[figLoc 'FGSC' saveTit ' PI Gating']);
    close all;
end
return;

function [errscore,flag,binsAboveThres,countsAboveThres] = otsufilter(intens,intensthres)
flag=0; % this flag will control warnings output to the command line

[counts,intensc] = hist(intens,255);
logcounts = log(max(1,counts));

isabovethres = intensc>intensthres & counts>10;
if sum(isabovethres)<4
    isabovethres = intensc>intensthres & counts>3; 
    flag=1;
end
    
if sum(isabovethres)<4
    flag=2;
end

binsAboveThres = intensc(isabovethres);
countsAboveThres = counts(isabovethres);

pfit=polyfit(binsAboveThres-intensthres,logcounts(isabovethres),1); % fit the counts data that is above the intensity threshold

predlogcounts = (pfit(2)+pfit(1)*(binsAboveThres-intensthres));
res = (predlogcounts-logcounts(isabovethres)).^2;

errscore = sum(res)/sum((abs(logcounts(isabovethres)-mean(logcounts(isabovethres)))).^2); % calculate error score for data above gate
return;