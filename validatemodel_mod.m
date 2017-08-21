function [Stats,T,hf] = validatemodel_mod(obj,varargin)
% VALIDATEMODEL Validates the quality of the credit scorecard model 
%
% Syntax:
%
%   Stats = validatemodel(sc)
%   Stats = validatemodel(sc,data)
%   [Stats,T] = validatemodel(sc,param1,value1,...)
%   [Stats,T,hf] = validatemodel(sc,param1,value1,...)
%
% Description:
%
%   VALIDATEMODEL Validates the quality of the credit scorecard model. The 
%   user can choose from a list of supported validation methods.
%
% Input Parameters 
%
%    sc  - Creditscorecard object
%
%    data - (Optional) Validation data, specified as a MATLAB table where
%       each row corresponds to individual observations. The data must 
%       contain columns for each of the predictors in the credit scorecard 
%       model. The columns of data can be any of the following: numeric,
%       logical, cell array of strings, character array, categorical. The
%       table must also contain a binary response variable.
%
% Optional Input Parameter Name/Value Pairs:
%
%   'AnalysisLevel' - String containing the type of analysis level to
%      perform when computing the validation statistics. Values are:
%
%              o 'Scores'   (default) Returns statistics (see output 'Stats')
%                           at the observation level. Scores are sorted from 
%                           riskiest to safest, and duplicates are removed.
%
%              o 'Deciles'  Returns the statistics (see output 'Stats') at 
%                           the decile level. Scores are sorted from riskiest 
%                           to safest, and duplicates are removed. Then, 
%                           observations are split into ten deciles (10%, 
%                           20%, ..., 100%) with corresponding statistics.                           
%
%   'Plot' - String or cell array containing the name of the plots to be
%      displayed. Possible values are:
%
%              o 'None'  (default) No plot is displayed.
%
%              o 'CAP'   Cumulative Accuracy Profile. Plots the fraction of 
%                        borrowers up to score "s" versus the fraction of 
%                        defaulters up to score "s" ('PctObs' versus 
%                        'Sensitivity' columns of T optional output argument).
%                        For more details, see Cumulative Accuracy Profile 
%                        (CAP) in the notes below.
%
%              o 'ROC'   Receiver Operating Characteristic. Plots the 
%                        fraction of non-defaulters up to score "s" versus 
%                        the fraction of defaulters up to score "s" 
%                        ('FalseAlarm' versus 'Sensitivity' columns of T 
%                        optional output argument). For more details, see 
%                        Receiver Operating Characteristic (ROC) in the 
%                        notes below.
%
%              o 'KS'    Kolmogorov-Smirnov statistic. Plots each score "s"
%                        versus the fraction of defaulters up to score "s",
%                        and also versus the fraction of non-defaulters up 
%                        to score "s," that is, 'Scores' versus both 
%                        'Sensitivity' and 'FalseAlarm' columns of the 
%                        optional output argument T. For more details, see 
%                        Kolmogorov-Smirnov statistic (KS) in the notes 
%                        below.
%      
%
% Output:
%
%  Stats  -  A 4-by-2 table containing the validation measures. The first 
%     column, 'Measure', contains the names of the following measures: 1) 
%     The Accuracy Ratio (AR), 2) the Area Under the ROC curve (AUROC),
%     3) the K-S statistic, and the 4) the K-S score.
%     The second column, 'Value', contains the values corresponding to these
%     measures.
%
%  T  - (Optional) Validation statistics data, returned as a N-by-9 summary 
%     table, sorted by score, from riskiest to safest. When 'AnalysisLevel'
%     is set to 'Deciles', N is equal to 10. Otherwise, N is equal to the 
%     total number of unique scores, i.e. scores without duplicates.
%     The table T contains the following information, in this order:
%
%         o 'Scores'       Scores sorted from riskiest to safest. The data 
%                          in this row corresponds to all observations up 
%                          to, and including, the score in this row.
%         o 'ProbDefault'  Probability of default for observations in this 
%                          row. For deciles, the average probability of 
%                          default for all observations in the given decile 
%                          is reported.
%         o 'TrueBads'     Cumulative number of "Bads" up to, and including,
%                          the corresponding score.
%         o 'FalseBads'    Cumulative number of "Goods" up to, and including, 
%                          the corresponding score.
%         o 'TrueGoods'    Cumulative number of "Goods" above the 
%                          corresponding score.
%         o 'FalseGoods'   Cumulative number of "Bads" above the 
%                          corresponding score.
%         o 'Sensitivity'  Fraction of defaulters (or cumulative number of 
%                          "bads" divided by total number of "bads"), this 
%                          is the distribution of "bads", up to, and 
%                          including, the corresponding score.
%         o 'FalseAlarm'   Fraction of non-defaulters (or cumulative number 
%                          of "goods" divided by total number of "goods"), 
%                          this is the distribution of "goods", up to, and 
%                          including, the corresponding score.
%         o 'PctObs'       Fraction of borrowers, or cumulative number of 
%                          observations, divided by total number of 
%                          observations, up to, and including, the 
%                          corresponding score.
%
%  hf  -  Handle or array of handles to the figure(s) of the measures being 
%     plotted. When 'Plot' is set to 'none', hf is an empty array.
%
% Notes:
%
% o The Cumulative Accuracy Profile (CAP), also known as the Gini curve, 
%   Power curve or Lorenz curve, the CAP is (generally) a concave curve. 
%   The scores of given observations are sorted from riskiest to safest. 
%   Then, for a given fraction M (0% to 100%) of the total borrowers, the 
%   height of the CAP curve is the fraction of defaulters whose scores are 
%   less than or equal to the maximum score of the fraction M, also known 
%   as 'Sensitivity'. The area under the CAP curve, known as the AUCAP, is 
%   then compared to that of the perfect or "ideal" model, leading to the 
%   definition of a summary index known as the Accuracy Ratio (AR) or the 
%   Gini Coefficient:
%   
%               AR = A_r/A_p,
%
%   Where A_r is the area between the CAP curve and the diagonal, and A_p 
%   is the area between the perfect model and the diagonal. The diagonal
%   represents a "random" model, one in which scores are assigned randomly 
%   and therefore the proportion of defaulters and non-defaulters is 
%   independent of the score. The perfect model is the model for which all 
%   defaulters are assigned the lowest scores, and therefore, perfectly 
%   discriminates between defaulters and non-defaulters. Thus, the closer 
%   to unity AR is, the better the scoring model.
%
% o The Receiver Operating Characteristic (ROC) is similar to the CAP. 
%   To find the ROC curve, we compute the proportion of defaulters up to a 
%   given score "s". This proportion is known as the True Positive Rate 
%   (TPR), or the "Sensitivity". Additionally, we also compute the proportion 
%   of non-defaulters up to score "s", also known as the False Positive Rate 
%   (FPR), or the 'False Alarm Rate'. The ROC curve is simply the plot of 
%   the 'Sensitivity' vs. the 'False Alarm Rate'. Computing the ROC curve is
%   similar to computing the equivalent of a confusion matrix at the 
%   score level.
%   Similar to the CAP, ROC has a summary index known as the area under the
%   ROC curve (AUROC). The closer to unity, the better the rating model. The
%   accuracy ratio (AR) is related to the area under the curve by the 
%   following formula:
%   
%               AR = 2*AUROC - 1.
%
% o The Kolmogorov-Smirnov (K-S) plot, also known as the fish-eye graph, 
%   is a common statistic used to measure the predictive power of 
%   scorecards. The KS plot shows the distribution of defaulters and the 
%   distribution of non-defaulters on the same plot. For the distribution 
%   of defaulters, each score "s" is plotted versus the proportion of 
%   defaulters up to "s", or "Sensitivity." For the distribution of 
%   non-defaulters, each score "s" is plotted versus the proportion of 
%   non-defaulters up to "s," or "False Alarm." The statistic of interest 
%   is called the KS statistic and is the maximum difference between these 
%   two distributions ("Sensitivity" minus "False Alarm"). The score at
%   which this maximum is attained is also of interest.
%
% References:
%
% o Basel Committee on Banking Supervision: Studies on the Validation of 
%   Internal Rating Systems, Working Paper No. 14, February 2005
%
% o Credit Risk Scorecards: Development and Implementation Using SAS, 
%   Refaat, M., Lulu.com, 2011
%
% o Credit Risk Modeling Using Excel and VBA, Loeffler, G. and Posch, P. N.,
%   Wiley Finance, 2007
%

% Copyright 2015 The MathWorks, Inc.

    parser = inputParser;
    
    parser.addOptional('data',obj.RawData,@istable);
    
    parser.addParameter('AnalysisLevel','Scores',@(x)ismember(lower(x),...
            {'deciles','scores'}));
        
    parser.addParameter('Plot','None',@(x)((ischar(x) || iscell(x)) ...
        && all(ismember(lower(x),{'cap','roc','ks','k-s','none'}))));
    
    try
        parser.parse(varargin{:});
    catch MException
        if strcmpi(MException.identifier,'MATLAB:InputParser:ArgumentFailedValidation')
            error(message('finance:creditscorecard:validatemodel:InvalidInputArguments'));
        else
            throw(MException);
        end
    end
    
    data = parser.Results.data;
    plotID = parser.Results.Plot;
    AnalysisLevel = parser.Results.AnalysisLevel;
    
    data = obj.PreprocessPredictorData(data);
    
    hasResponse = ismember(obj.ResponseVar,data.Properties.VariableNames);
    if ~hasResponse
       error(message('finance:creditscorecard:validatemodel:MissingResponse'));
    end

    % Remove rows with missing response
    ResponseMissingDataInd = isnan(double(data.(obj.ResponseVar)));
    data(ResponseMissingDataInd,:) = [];
    if height(data)==0
       error(message('finance:creditscorecard:validatemodel:ResponseIsMissingAll'));
    end
    
    Scores      = score_mod(obj,data);
    ProbDefault = probdefault_mod(obj,data);
    RespVar     = data.(obj.ResponseVar);
    RespValues  = validateResponseValues(obj,obj.RawData.(obj.ResponseVar));
    [~,RespVarMap] = MapResponseVar(obj,RespVar,RespValues);
    NumDef      = sum(RespVarMap == 0);
    NumObs      = numel(Scores);

    % Get first score data for KS plot
    [~,iWorst]  = max(ProbDefault);
    WorstScore  = Scores(iWorst);
    WorstData   = RespVarMap(Scores == WorstScore);
    WorstGoods  = sum(WorstData);
    WorstBads   = length(WorstData) - WorstGoods;
    
    % Compute table of data or decile level statistics
    T = getValidationTable(ProbDefault,Scores,...
        RespVarMap,AnalysisLevel,NumDef,NumObs);

    % Compute relevant measures
    Stats = getValidationStatistics(T,NumDef,NumObs);

    % Show plots
    hf = plotValidationCurves(plotID,T,Stats,NumDef,NumObs,...
       WorstScore,WorstBads,WorstGoods);

end

%% Helper functions

function T = getValidationTable(ProbDefault,Scores,RespVarMap,AnalysisLevel,NumDef,NumObs)
    % This helper function computes the statistics used to plot the curves.
    % The reported statistics are either at the score or the decile level.
    % These are the Scores, sorted from riskiest to safest, the 
    % corresponding Probabilities of Default, True Bads, False Bads, True 
    % Goods, False Goods, Sensitivity, False Alarm Rate and Percentage 
    % Observations.

    if strcmpi(AnalysisLevel,'deciles')
       
        % Create a Numeric container for ProbDefault against mapped
        % response. Apply Equal Frequency with 10 bins.
        % The container sorts from low probabilities to high probabilities,
        % and the cut points use the "left-edge" convention. This works
        % well, because we need a "right-edge" convention sorting from
        % high to low probabilities. So, after binning, we simply need to
        % flip the order of the cut points and the order of the frequency
        % table.
        dbc = internal.finance.binning.container.Numeric(ProbDefault,...
           RespVarMap,'ResponseOrder',[0 1]);
        % Create an equal frequency algorithm (try to) find the 10 equal
        % frequency bins
        ef = internal.finance.binning.algorithm.EqualFrequency;
        ef = ef.setOptions('NumBins',10);
        dbc = ef.runAlgorithm(dbc);

        % Get cumulative frequency table; note we flip, as per discussion
        % above.
        DecileCumFTPD = cumsum(flipud(dbc.getFrequencyTable));

        % We need to map cut points from pd's to scores, including the best
        % possible score, for the table cut points.
       
        % Find cut points in the order we want them (flipped) and including
        % the min pd; keep location of min pd (best score)
        [MinPD,MinPDInd] = min(ProbDefault);
        DecileCPPD = [flipud(dbc.getCutPoints);MinPD];

        % The number of bins after binning might not be 10, if there are too
        % few possible scores, or if only a few scores are assigned to most
        % of the data.
        NumDecileBinsPD = length(DecileCPPD);
       
        % Currently, the cut points from the container are not
        % full-precision, so we need to find the full-precision pd's before
        % mapping and before aggregating pd's.
        DecileCPPDFullPrecision = zeros(size(DecileCPPD));
        DecileCPScore = zeros(size(DecileCPPD));
        for ii=1:NumDecileBinsPD-1
           [~,ind] = min(abs(DecileCPPD(ii)-ProbDefault));
           DecileCPPDFullPrecision(ii) = ProbDefault(ind);
           DecileCPScore(ii) = Scores(ind);
        end
        % Last DecileCPPD is full-precision, by construction
        DecileCPPDFullPrecision(end) = DecileCPPD(end);
        DecileCPScore(end) = Scores(MinPDInd);

        % Compute pd's for each bin as the average pd of the data points in
        % the bin. Note that this implementation is actually a weighted
        % average, since the average is computed over all observations in
        % the bin.
        DecilePD = zeros(NumDecileBinsPD,1);
        DecilePD(1) = mean(ProbDefault(ProbDefault>=DecileCPPDFullPrecision(1)));
        for ii=2:NumDecileBinsPD
           DecilePD(ii) = mean(ProbDefault(DecileCPPDFullPrecision(ii-1)>ProbDefault &...
              ProbDefault>=DecileCPPDFullPrecision(ii)));
        end

        % Set variables used to finish the validation table below, outside
        % the decile v. score branching
        sScores = DecileCPScore;
        sProbDefault = DecilePD;
        TrueBads  = DecileCumFTPD(:,1);
        FalseBads = DecileCumFTPD(:,2);
        CumObs = TrueBads + FalseBads;
        PctObs = CumObs/CumObs(end);

    else
        
        [sProbDefault,Ind] = sort(ProbDefault,'descend');
        sScores = Scores(Ind);
        sRespVarMap = RespVarMap(Ind);

        TrueBads  = cumsum(~sRespVarMap); % CumBads
        FalseBads = cumsum(sRespVarMap); % CumGoods
        CumObs    = TrueBads + FalseBads;
        PctObs    = CumObs./NumObs;

        % Remove duplicates from sorted scores and statistics
        uInd = [sScores(1:end-1) ~= sScores(2:end);true];
        sScores = sScores(uInd);
        sProbDefault = sProbDefault(uInd);
        TrueBads  = TrueBads(uInd);
        FalseBads = FalseBads(uInd);
        PctObs    = PctObs(uInd);
        
    end
    
    % Compute the True Goods, False Goods, Sensitivity and False Alarm Rate
    TrueGoods   = FalseBads(end) - FalseBads;
    FalseGoods  = TrueBads(end) - TrueBads;
    Sensitivity = TrueBads./NumDef;
    FalseAlarm  = FalseBads./(NumObs-NumDef);
    
    T = table(sScores,sProbDefault,TrueBads,FalseBads,TrueGoods,FalseGoods,...
        Sensitivity,FalseAlarm,PctObs,'VariableNames',{'Scores','ProbDefault',...
        'TrueBads','FalseBads','TrueGoods','FalseGoods','Sensitivity',...
        'FalseAlarm','PctObs'});

end

function Stats = getValidationStatistics(T,NumDef,NumObs)

    % CAP
    x = [0;T.PctObs];
    y = [0;T.Sensitivity];
    
    aucap = 0.5*(y(1:end-1)+y(2:end))'*diff(x);
    ar = (aucap-0.5)/((1-NumDef/NumObs/2)-0.5);

    % ROC
    x = [0;T.FalseAlarm];
    % same y as CAP
    
    auroc = 0.5*(y(1:end-1)+y(2:end))'*diff(x);
    
    % KS
    [KSValue,KSInd] = max(T.Sensitivity-T.FalseAlarm);
    KSScore = T.Scores(KSInd);
    
    % Create Stats table output
    Measure = {'Accuracy Ratio','Area under ROC curve','KS statistic',...
            'KS score'}';
    Value  = [ar;auroc;KSValue;KSScore];
    
    Stats = table(Measure,Value);

end

%% Plot functions

function hf = plotValidationCurves(plotID,T,Stats,NumDef,NumObs,...
   WorstScore,WorstBads,WorstGoods)
    
    xPerfMdl = [0; NumDef/NumObs; 1];
    yPerfMdl = [0; 1; 1];
    Color1   = [1 1 0.6]; 
    Color2   = [0.8 1 1];
    
    % CAP data
    PctObs = [0;T.PctObs];
    TrueBadsRate = [0;T.Sensitivity];
    
    % ROC data
    FalseBadsRate = [0;T.FalseAlarm];
    
    % KS data
    Scores   = T.Scores;
    CumBads  = T.Sensitivity;
    CumGoods = T.FalseAlarm;
    
    % Additional KS adjustments and information
    
    % Adjust first point for deciles. Strictly speaking, the adjustment is
    % for any situation where the first value in Scores is not the minimum
    % observed score (which is the usual situation with deciles), in which
    % case we want to add data for the left-most values of the KS plot.
    if Scores(1) ~= WorstScore
       Scores   = [WorstScore;Scores];
       CumBads  = [WorstBads/NumDef; CumBads];
       CumGoods = [WorstGoods/(NumObs-NumDef); CumGoods];
    end
    
    % Get values for text in CAP plot
    [~,IndAR] = ismember('Accuracy Ratio',Stats.Measure);
    ar = Stats.Value(IndAR);
    
    % Get values for text in ROC plot
    [~,IndAUROC] = ismember('Area under ROC curve',Stats.Measure);
    auroc = Stats.Value(IndAUROC);
    
    % Get values for text in KS plot
    [~,IndKSValue] = ismember('KS statistic',Stats.Measure);
    KSValue = Stats.Value(IndKSValue);
    [~,IndKSScore] = ismember('KS score',Stats.Measure);
    KSScore = Stats.Value(IndKSScore);
    KSInd = find(KSScore==Scores,1,'first');
    
    if ischar(plotID)
        plotID = {plotID};
    end
        
    plotID = unique(plotID,'stable');
    
    hFig = cell(numel(plotID),1);
    
    for i = 1 : numel(plotID)
        % CAP
        if strcmpi(plotID{i},'cap')
            hFig{i} = plotCAP(PctObs,TrueBadsRate,xPerfMdl,yPerfMdl,...
                Color1,Color2,ar);
        end
        
        % ROC
        if strcmpi(plotID{i},'roc')
            hFig{i} = plotROC(FalseBadsRate,TrueBadsRate,Color1,auroc);
        end
        
        % KS
        if any(strcmpi(plotID{i},{'ks','k-s'}))
            hFig{i} = plotKS(Scores,CumBads,CumGoods,KSScore,KSValue,KSInd);
        end
    end
    
    hf = [hFig{:}];
end


function hcap = plotCAP(x,y,xPerfMdl,yPerfMdl,Color1,Color2,Measure)
    hcap = figure;
    hax  = axes('Parent',hcap);

    xLimits = get(hax,'XLim');
    yLimits = get(hax,'YLim');
    

    hp1 = fill([x;flipud(x)],[x;flipud(y)],Color1);
    set(hp1,'Parent',hax,'EdgeColor','none')
    hold(hax,'on')
    hp2 = fill([x;flipud(xPerfMdl)],[y;flipud(yPerfMdl)],Color2);
    set(hp2,'Parent',hax,'EdgeColor','none')

    box(hax,'on')
    plot(hax,x,y,'k-')
    plot(hax,x,x,'k--')
    plot(hax,xPerfMdl,yPerfMdl,'k')
    TextInPlot = sprintf('AR = %0.3f',Measure);
    text(0.6*diff(xLimits),0.2*diff(yLimits),TextInPlot,...
        'HorizontalAlignment','Center','Parent',hax)
    xlabel('Fraction of borrowers')
    ylabel('Fraction of defaulters')
    title(hax,'Cumulative Accuracy Profile (CAP) curve')

end


function hroc = plotROC(x,y,Color,Measure)
    hroc = figure;
    hax  = axes('Parent',hroc);

    xLimits = get(hax,'XLim');
    yLimits = get(hax,'YLim');
    X = [x;1;1;0];
    Y = [y;1;0;0];
    
    hp1 = fill(X,Y,Color);
    set(hp1,'Parent',hax,'EdgeColor','none')
    hold(hax,'on')
    plot(hax,x,y,'k-')
    box(hax,'on')
    TextInPlot = sprintf('AUROC = %0.3f',Measure);
    text(0.6*diff(xLimits),0.2*diff(yLimits),TextInPlot,...
        'HorizontalAlignment','Center','Parent',hax)
    xlabel('Fraction of non-defaulters')
    ylabel('Fraction of defaulters')
    title(hax,'Receiver Operating Characteristic (ROC) curve')
    
end


function hks = plotKS(Scores,CumBads,CumGoods,KSScore,KSValue,KSInd)
    hks = figure;
    hax = axes('Parent',hks);
    
    plot(hax,Scores,CumBads,Scores,CumGoods);
    if all(diff(Scores) < 0)
        set(hax,'XDir','Reverse');
    end
    hl = legend({'Cumulative Bads','Cumulative Goods'},'location','Best',...
        'AutoUpdate','Off');
    set(hl,'Parent',hks)
    xlabel('Score (Riskiest to Safest)')
    ylabel('Cumulative Probability')

    hold(hax,'on')
    xLimits = get(hax,'XLim');
    yLimits = get(hax,'YLim');

    plot(hax,[KSScore KSScore],yLimits,'k:')
    plot(hax,[xLimits(1) KSScore],[CumBads(KSInd) CumBads(KSInd)],'k:')
    plot(hax,[xLimits(1) KSScore],[CumGoods(KSInd) CumGoods(KSInd)],'k:')

    TextInPlot = sprintf('K-S %3.1f%%, at %g',KSValue*100,KSScore);
    text((xLimits(1)+KSScore)/2,(CumGoods(KSInd)+CumBads(KSInd))/2,TextInPlot,...
        'HorizontalAlignment','Center','Parent',hax)
    
    title(hax,'K-S Plot')
    hold(hax,'off')
end


