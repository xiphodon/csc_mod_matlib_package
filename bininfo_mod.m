function [bi,bm,mv] = bininfo_mod(obj,varargin)
% BININFO Return predictor's bin information.
%
% Syntax:
%
%   bi = bininfo(sc,PredictorName)
%   bi = bininfo(sc,PredictorName,param1,val1,...)
%   [bi,bm] = bininfo(sc,PredictorName,param1,val1,...)
%   [bi,bm,mv] = bininfo(sc,PredictorName,param1,val1,...)
%
% Description:
%
%   For the predictor specified in PredictorName, BININFO provides
%   information at bin level, such as frequencies of Goods and Bads, Odds,
%   Weight Of Evidence (WOE), etc. It also returns the binning map or rules
%   in the form of a vector of cut points for numeric predictors, or a
%   table of category groupings for categorical predictors.
%
% Input Arguments:
%
%   sc - creditscorecard object.
%
%   PredictorName - String containing the predictor's name.
%
% Optional Input Parameter Name/Value Pairs:
%
%   'Statistics'   List of statistics to include in the bin information.
%                  By default, 'Odds','WOE','InfoValue' are reported. 
%                  'Entropy', when requested, is also reported.
%                  The input must be a cell containing the name(s) of 
%                  the desired statistics:
%
%                   o 'Odds'       Ratio of Goods over Bads.
%
%                   o 'WOE'        Statistic used to measure the deviation
%                                  between the distribution of Goods and 
%                                  Bads.
%
%                   o 'InfoValue'  Closely tied to the WOE, it is a 
%                                  statistic used to determine how strong a 
%                                  predictor is to use in the fitting
%                                  model. It measures how strong the 
%                                  deviation between the distributions of 
%                                  Goods and Bads.
%
%                   o 'Entropy'    Measure of unpredictability contained in
%                                  the bins. Used as a technique to 
%                                  validate a risk model.
%
%
%   'Totals'       'On' (default) | 'Off'. Include a row of totals at the 
%                  bottom of the information table.
%
% Output:
%
%   bi - Bin information table. It contains one row per bin and a row of
%      totals. The columns contain bin descriptions, frequencies of goods
%      and bads, and bin statistics.
%
%   bm - Binning map or rules in the form of a vector of cut points for
%      numeric predictors, or a table of category groupings for
%      categorical predictors.
%
%   mv - Numeric array containing the minimum and maximum values, as set 
%       (or defined) by the user. This output argument is set to an empty 
%       array for categorical predictors.
%
% Algorithms:
%
%   Suppose the predictor's data takes on M possible values b1, ..., bM.
%   For binned data, M is a small number. The response takes on two values,
%   'Good' and 'Bad'. The frequency table of the data is given by
%
%              Good     Bad    Total
%
%         b1:   n11     n12     n1
%         b2:   n21     n22     n2
%                 ...
%         bM:   nM1     nM2     nM
%
%      Total:  nGood   nBad   nTotal
%
%   The Weight Of Evidence (WOE) is defined for each data value bi as
%
%         WOE(i) = log((ni1/nGood)/(ni2/nBad)).
%
%   If we define
%
%         pGood(i) = ni1/nGood, pBad(i) = ni2/nBad
%
%   then pGood(i) is the proportion of 'good' observations that take on the
%   value bi, and similarly for pBad(i). In other words, pGood(i) gives the
%   distribution of good observations over the M observed values of the
%   predictor, and similarly for pBad(i). With this, an equivalent formula
%   for the WOE is
%
%         WOE(i) = log(pGood(i)/pBad(i)).
%
%   Using the same frequency table, the Odds for row i are defined as
%
%         Odds(i) = ni1 / ni2,
%
%   and the Odds for the sample are defined as
%
%         OddsTotal = nGood / nBad.
%
%   For each row i we can also compute its contribution to the total
%   Information Value, given by
%
%         InfoValue(i) = (pGood(i) - pBad(i)) * WOE(i),
%
%   and the total Information Value is simply the sum of all the InfoValue(i)
%   terms.
%
%  Likewise, for each row i, we can compute its contribution to the total
%  Entropy, given by
%
%         Entropy(i) = -1/log(2)*(ni1/ni*log(ni1/ni)+ni2/ni*log(ni2/ni),
%
%  and the total Entropy is simply the weighted sum of the row entropies,
%
%         Entropy = sum(ni/nTotal * Entropy(i)), i = 1...M
%
%
% Notes:
%
%   o The Weight Of Evidence (WOE) is only defined at a bin level, and
%     there is no definition for a 'Total WOE'. The bin information table
%     sets the Total WOE to NaN.
%
%   o It is recommended to avoid having bins with frequencies of zero,
%     because they lead to infinite or undefined (NaN) statistics. Use
%     MODIFYBINS or AUTOBINNING to modify bins.
%
%   o The computation of Information Value does ignore bins that are empty
%     (neither good nor bad observations), because it uses a 'nansum' to
%     compute the total Information Value. However, bins with only good or
%     only bad observations do lead to an infinite Information Value.
%     Consider modifying the bins in those cases (see previous point).
%
%  o The values in the third output argument are the minimum and maximum 
%    allowed values for the predictor, set by default to -Inf and Inf. The
%    user can modify these using the optional inputs "MinValue" and 
%    "MaxValue" of the method MODIFYBINS. These minimum and maximum values 
%    are the same as the values in the bin labels in the first output of 
%    BININFO. Specifically, the left edge of the first bin, and the right 
%    edge of the last bin. To see the minimum and maximum predictor values 
%    observed in the data, use the second output of PREDICTORINFO.
%
% See also MODIFYBINS, AUTOBINNING, PLOTBINS, BINDATA.

% Copyright 2014 - 2016 The MathWorks, Inc.
    data   = obj.RawData;        
    parser = inputParser;

    defaultStats = {'odds','woe','infoval'};

    parser.addRequired('PredictorVar',@(x)(ischar(x) && size(x,1) == 1));
    parser.addParameter('Statistics',defaultStats,@(x)(ischar(x)||iscell(x)));
    parser.addParameter('Totals','On',@(x)ismember(lower(x),{'on','off'}));                

    if nargin == 1
        error(message('finance:creditscorecard:bininfo:NotEnoughInputArguments'))
    else
        try
            parser.parse(varargin{:});
        catch MException
            if ~ismember(varargin{1},setdiff(obj.PredictorVars,...
                        {obj.IDVar obj.ResponseVar}))
                    
                    error(message('finance:creditscorecard:bininfo:MissingPredictorName'));
                    
            elseif strcmpi(MException.identifier,...
                        'MATLAB:InputParser:ParamMissingValue')
                    
                    error(message('finance:creditscorecard:bininfo:IncompatibleNameValuePair'));
            else
                throw(MException)
            end
        end
    end

    PredictorVar = parser.Results.PredictorVar;
    StatsList    = parser.Results.Statistics;
    Totals       = parser.Results.Totals;
    
    try 
        PredictorVar = obj.PreprocessPredictorVar(PredictorVar);
        
    catch MException
        % Case where the inputParser passes, but both the PredictorVar 
        % and the value of the parameter is missing: 
        % e.g. bininfo(sc,'Statistics')
        throw(MException);
    end

    [obj,Ind] = obj.getContainerIndex(PredictorVar{:});
    BinLabels = obj.BinContainers{Ind}.getLabels;
    [FreqTab, NanFreqTab] = obj.BinContainers{Ind}.getFrequencyTable();

    if ischar(StatsList)
        StatsList = {StatsList};
    end
    
    [okIV,idxInfoVal] = ismember('infovalue',lower(StatsList));
    [okEnt,idxEnt]  = ismember('entropy',lower(StatsList));
    
    if okIV
        StatsList{idxInfoVal} = 'InfoVal';
    end
    
    if okEnt
        StatsList{idxEnt} = 'Entropy';
    end
    
    Stats = obj.BinContainers{Ind}.getStatistics('StatsList',StatsList);   

    [StatNames,StatValues] = formatStats(Stats,StatsList);

    if any(NanFreqTab ~= 0)
        bi = table(BinLabels,[FreqTab(:,1);NanFreqTab(1)],[FreqTab(:,2);NanFreqTab(2)],StatValues{:,:},...
            'VariableNames',['Bin' 'Good' 'Bad' StatNames(:)']);      
    else
        bi = table(BinLabels,FreqTab(:,1),FreqTab(:,2),StatValues{:,:},...
            'VariableNames',['Bin' 'Good' 'Bad' StatNames(:)']);   
    end
        
    if strcmpi(Totals,'on')
        C = cellfun(@(c)isnumeric(bi.(c)),bi.Properties.VariableNames,...
            'UniformOutput',false);
        C = [C{:}];

        TotTable  = array2table(zeros(size(C)),'VariableNames',...
                    bi.Properties.VariableNames);

        TotTable.Bin = {'Totals'}; 
        bi(end+1,:)    = TotTable;

        for i = 1 : numel(C)
            if C(i)
                if strcmpi(bi.Properties.VariableNames{i},'WOE')
                    bi(end,i).(bi.Properties.VariableNames{i}) = NaN;

                elseif strcmpi(bi.Properties.VariableNames{i},'Odds')
                    bi(end,i).(bi.Properties.VariableNames{i}) = ...
                        bi.Good(end)/bi.Bad(end);
                    
                elseif strcmpi(bi.Properties.VariableNames{i},'InfoValue')
                    bi(end,i).(bi.Properties.VariableNames{i}) = Stats.infoval;
                    
                elseif strcmpi(bi.Properties.VariableNames{i},'Entropy')
                    bi(end,i).(bi.Properties.VariableNames{i}) = Stats.entropy;    
                    
                else 
                    bi(end,i).(bi.Properties.VariableNames{i}) = ...
                     nansum(bi.(bi.Properties.VariableNames{i}));
                 
                end
            end
        end
    end

    % Get cut points or category groupings
    if isnumeric(data.(PredictorVar{:}))
        bm = obj.BinContainers{Ind}.getCutPoints();
    else
        bm = obj.BinContainers{Ind}.getCatGrouping();
    end
    
    % If the predictor is numeric and the user defined Min and Max valued
    % are queried, return them as third output
    if isnumeric(data.(PredictorVar{:})) && (nargout == 3)
        mv = [obj.BinContainers{Ind}.MinValue obj.BinContainers{Ind}.MaxValue];
    else
        mv = [];
    end


    function [Snames,Svalues] = formatStats(S,Slist)
        % Returns formatted statistics names and values, ready to use
        % to create the bin information table.
        names   = fieldnames(S);
        values  = cellfun(@(c)S.(c),names,'UniformOutput',false);
        
        Sfields   = fieldnames(S);
        IVnames   = {'infoval' 'infovalbins'};
        Entnames  = {'entropy' 'entropybins'};
        [~,Indx]  = ismember(lower(Slist),Sfields);
        [~,iiv1]  = ismember(IVnames,lower(Slist));
        [okIv,iiv2] = ismember(IVnames,Sfields);
        [~,ient1]   = ismember(Entnames,lower(Slist));
        [okE,ient2] = ismember(Entnames,Sfields);
        
        % If 'infoval' is a field in 'S', swap its index with that of 
        % 'infovalbins', since these latter (bin-level values) will be used 
        % in the table. The same applies for 'entropybins'.
        % The same applies to 'entropybins'.
        idxIV = strcmpi('infoval',IVnames);
        idxE  = strcmpi('entropy',Entnames);
        
        if okIv(idxIV) 
            Indx(iiv1(idxIV)) = iiv2(strcmpi('infovalbins',IVnames));
        end
        
        if okE(idxE) 
            Indx(ient1(idxE)) = ient2(strcmpi('entropybins',Entnames));
        end
        
        % Keep the same mapping between the names in 'Snames' and the desired
        % statistics values in 'Svalues'
        Snames  = Slist; 
        Svalues = values(Indx);
        
        if ischar(Snames)
            Snames = {Snames};
        end
              
        % Rename the InfoVal and Entropy columns, and format names 
        Snames(strcmpi(Snames,'odds')) = {'Odds'};
        Snames(strcmpi(Snames,'infoval')) = {'InfoValue'};
        Snames(strcmpi(Snames,'entropy')) = {'Entropy'};
        Snames(strcmpi(Snames,'woe')) = {'WOE'};
    end

end
    
    