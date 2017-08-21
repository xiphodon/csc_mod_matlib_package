function obj = autobinning_mod(obj,varargin)
% AUTOBINNING Perform automatic binning of given predictors.
%
% Syntax:
%
%   sc = autobinning(sc)
%   sc = autobinning(sc,PredictorNames)
%   sc = autobinning(sc,PredictorNames,param1,val1,...)
%
% Description:
%
%   AUTOBINNING performs automatic binning of the predictors given in
%   PredictorNames. See the Algorithms section for a description of the
%   supported binning algorithms.
%
%   Automatic binning finds binning maps or rules to bin numeric data or to
%   group categories of categorical data. The binning rules are stored in
%   the creditscorecard object. To apply the binning rules to the 
%   creditscorecard data, or to a new dataset, use the method BINDATA .
%
% Input Arguments:
%
%   sc - creditscorecard object.
%
%   PredictorNames -  A string or a cell array of strings containing the
%      name of the predictor or predictors that need to be automatically
%      binned. When no predictor name is passed, all predictors in the
%      'PredictorVars' property of the creditscorecard object are binned.
%
% Optional Input Parameter Name/Value Pairs:
%
%   'Algorithm'    A string indicating the algorithm to be used (same
%                  algorithm for all predictors in PredictorNames). 
%                  Possible values are:
%
%                  o 'Monotone'         (default) Monotone Adjacent Pooling
%                                       Algorithm (MAPA), also know as 
%                                       Maximum Likelihood Coarse 
%                                       Classifier. Supervised algorithm 
%                                       that aims to find bins with a 
%                                       monotone Weight-Of-Evidence
%                                       (WOE) trend.
%                                       This algorithm assumes that only 
%                                       neighboring attributes can be 
%                                       grouped. Thus, for categorical 
%                                       predictors, categories are 
%                                       preliminarily sorted before applying 
%                                       the algorithm.
%
%                  o 'EqualFrequency'   Unsupervised algorithm that divides 
%                                       the data into a predetermined number 
%                                       of bins that contain approximately 
%                                       the same number of observations. 
%                                       This algorithm is also known as 
%                                       "equal height" or "equal depth".
%                                       For categorical predictors,
%                                       categories are sorted before
%                                       applying the algorithm.
%
%                  o 'EqualWidth'       Unsupervised algorithm that divides 
%                                       the range of values in the domain
%                                       of the predictor variable into a 
%                                       predetermined number of bins of 
%                                       "equal width." For numeric data, 
%                                       the width is measured as the 
%                                       distance between bin edges. 
%                                       For categorical data, width is 
%                                       measured as the number of 
%                                       categories within a bin. For 
%                                       categorical predictors, categories 
%                                       are sorted before applying the 
%                                       algorithm.
%
%
%   'AlgorithmOptions' A cell array containing name-value pairs with
%                      algorithm options for the selected algorithm. For
%                      example, for the Equal Frequency algorithm, to set
%                      the number of bins to 7, use
%                      'AlgorithmOptions',{'NumBins',7}. Options are:
%
%                      For 'Monotone':
%                      o {'InitialNumBins',n}    Initial number (n) of bins 
%                                   (default is 10). n must be an integer 
%                                   number greater than 2. This option is 
%                                   used for numeric predictors only.
%
%                      o {'Trend','TrendOption'} Determines whether the 
%                                   Weight Of Evidence (WOE) monotonic 
%                                   trend is expected to be increasing or 
%                                   decreasing. The values for 
%                                   'TrendOptions' are:
%
%                                   'Auto'        (default) Automatically
%                                                 determines if the trend is 
%                                                 increasing or decreasing.
%
%                                   'Increasing'  Look for an increasing 
%                                                 WOE trend.
%
%                                   'Decreasing'  Look for a decreasing 
%                                                 WOE trend.
%                                   
%                                   Note: The value of the input 'Trend' 
%                                   does not necessarily reflect that of 
%                                   the resulting WOE curve. The parameter
%                                   'Trend' rather tells the algorithm to
%                                   "look for" an increasing or decreasing
%                                   trend, but the outcome may not show the
%                                   desired trend. For more information,
%                                   see the description of 'Monotone'
%                                   below, under 'Algorithms'.
%
%                      o {'SortCategories','SortOption'} For categorical 
%                                   predictors only. Used to determine how 
%                                   the predictor categories are sorted as 
%                                   preprocessing step before applying the 
%                                   algorithm. The values of 'SortOption'
%                                   are:
%
%                                   'Odds'   (default) The categories are 
%                                            sorted by order of increasing 
%                                            values of odds, defined as 
%                                            the ratio of "Good" to "Bad"
%                                            observations, for the given
%                                            category.
%                                   'Goods'  The categories are sorted by 
%                                            order of increasing values 
%                                            of "Good" count.
%                                   'Bads'   The categories are sorted by 
%                                            order of increasing values
%                                            of "Bad" count.
%                                   'Totals' The categories are sorted by 
%                                            order of increasing values  
%                                            of total number of observations
%                                            ("Good" plus "Bad").
%                                   'None'   No sorting is applied. 
%                                            The existing order of the 
%                                            categories is unchanged before 
%                                            applying the algorithm. (The 
%                                            existing order of the categories 
%                                            can be seen in the category 
%                                            grouping optional output from 
%                                            BININFO.)
%
%                                   For more information, see section "Sort 
%                                   Categories" in the documentation.
%
%                      For 'EqualFrequency':
%                      o {'NumBins',n}       Specifies the desired number 
%                                            (n) of bins (default is 5). n 
%                                            must be a positive number.
%
%                      o {'SortCategories','SortOption'}  See above, in 
%                                            Monotone.
%
%                                   Note: The existing order of the categories 
%                                   is unchanged before applying the 
%                                   algorithm. (The existing order of the 
%                                   categories can be seen in the category 
%                                   grouping optional output from BININFO.)
%                                   For more information, see section "Sort 
%                                   Categories" in the documentation.
%
%                      For 'EqualWidth':
%                      o {'NumBins',n}       Specifies the desired number 
%                                            (n) of bins (default is 5). n 
%                                            must be a positive number.
%
%                      o {'SortCategories','SortOption'}  See above, in 
%                                            Monotone.
%
%                                   Note: The existing order of the categories 
%                                   is unchanged before applying the 
%                                   algorithm. (The existing order of the 
%                                   categories can be seen in the category 
%                                   grouping optional output from BININFO.)
%                                   For more information, see section "Sort 
%                                   Categories" in the documentation.
%
%
%   'Display'   'On' | 'Off' (default). Indicator to display the information 
%               on status of the binning process at command line, specified 
%               using a string with value 'On' or 'Off'.
%
%
% Output:
%
%   sc - Updated creditscorecard object, containing the automatically
%      determined binning maps or rules (cut points or category groupings)
%      for the predictor(s).
%
% Algorithms:
%
%   Here is a description of the binning algorithms.
%
%   o Monotone: The 'Monotone' algorithm is an implementation of the 
%     Monotone Adjacent Pooling Algorithm (MAPA), also known as the Maximum
%     Likelihood Monotone Coarse Classifier (see Anderson or Thomas in the
%     references). 
%       
%     Preprocessing:
%
%     During the preprocessing phase, preprocessing of numeric predictors
%     consists in applying equal frequency binning, with the number of 
%     bins determined by the 'InitialNumBins' parameter (the default is 10 
%     bins). The preprocessing of categorical predictors consists in 
%     sorting the categories according to the 'SortCategories' criterion 
%     (the default is to sort by odds in increasing order). Sorting is 
%     not applied to ordinal predictors. See the definitions of "Sort 
%     Categories" or the description of AlgorithmOptions option for 
%     'SortCategories' in the documentation.
%
%     Main Algorithm:
%
%     The following example illustrates how the 'Monotone' algorithm arrives
%     at cut points for numeric data.
% 
%     Bin	        Good	Bad	  Iteration1	Iteration2	Iteration3	Iteration4
%    '[-Inf,33000)'  127	107	  0.543	 	 	 
%    '[33000,38000)' 194  	 90	  0.620	        0.683	 	 
%    '[38000,42000)' 135	 78   0.624	        0.662	 	 
%    '[42000,47000)' 164	 66	  0.645	        0.678	    0.713	 
%    '[47000,Inf]'   183	 56	  0.669	        0.700	    0.740	    0.766
% 
%    Initially, the numeric data is preprocessed with an equal frequency 
%    binning. In this example, for simplicity, only the five initial bins 
%    are used. The first column indicates the equal frequency bin ranges, 
%    and the second and third columns have the "Good" and "Bad" counts per 
%    bin.
%
%    Monotone finds break points based on the cumulative proportion of 
%    "Good" observations. In the 'Iteration1' column, the first value (0.543)
%    is the number of "Good" observations in the first bin (127), divided 
%    by the total number of observations in the bin (127+107). The second 
%    value (0.620) is the number of "Good" observations in bins 1 and 2, 
%    divided by the total number of observations in bins 1 and 2. And so 
%    forth. The first cut point is set where the minimum of this cumulative 
%    ratio is found, which is in the first bin in this example. This is the 
%    end of iteration 1.
% 
%    Starting from the second bin, cumulative proportions of "Good" 
%    observations are computed again. The second cut point is set where the 
%    minimum of this cumulative ratio is found. In this case, it happens to 
%    be in bin number 3, therefore bins 2 and 3 are merged.
% 
%    The algorithm proceeds the same way for two more iterations. In this 
%    particular example, in the end it only merges bins 2 and 3. The final 
%    binning has four bins with cut points at 33,000, 42,000, and 47,000.
% 
%    For categorical data, the only difference is that the preprocessing 
%    step consists in reordering the categories. Consider the following 
%    categorical data:
%
%    Bin	      Good	   Bad	   Odds
%    'Home Owner'  365	   177	  2.062
%    'Tenant'      307	   167    1.838
%    'Other'       131	    53	  2.474
% 
%    The preprocessing step, by default, sorts the categories by 'Odds'. 
%    (See the Sort Categories definition or the description of 
%    AlgorithmOptions option for 'SortCategories' for more information.) 
%    Then, it applies the same steps described above, shown in the following 
%    table:
% 
%    Bin	       Good   Bad	 Odds	Iteration1	Iteration2	Iteration3
%    'Tenant'	    307	  167	1.838	     0.648	 	 
%    'Home Owner'	365	  177	2.062	     0.661	     0.673	 
%    'Other'	    131	   53	2.472	     0.669	     0.683	     0.712
% 
%    In this case, the Monotone algorithm would not merge any categories. 
%    The only difference, compared with the data before the application of 
%    the algorithm, is that the categories are now sorted by 'Odds'.
% 
%    In both the numeric and categorical examples above, the implicit 
%    'Trend' choice is 'Increasing'. (See the description of AlgorithmOptions 
%    option for the 'Monotone' 'Trend' option.) If you set the trend to 
%    'Decreasing', the algorithm looks for the maximum (instead of the 
%    minimum) cumulative ratios to determine the cut points. In that case, 
%    at iteration 1, the maximum would be in the last bin, which would imply 
%    that all bins should be merged into a single bin. Binning into a single 
%    bin is a total loss of information and has no practical use. Therefore,
%    when the chosen trend leads to a single bin, the Monotone implementation 
%    rejects it, and the algorithm returns the bins found after the 
%    preprocessing step. This state is the initial equal frequency binning 
%    for numeric data and the sorted categories for categorical data. 
%    The implementation of the Monotone algorithm by default uses a heuristic 
%    to identify the trend ('Auto' option for 'Trend').
%
%
%   o Equal Frequency: Unsupervised algorithm that divides the data into a 
%     predetermined number of bins that contain approximately the same 
%     number of observations. 
%     Let v[1], v[2],..., v[N] be the sorted list of different values or 
%     categories observed in the data. Let f[i] be the frequency of v[i]. 
%     Let F[k] = f[1]+...+f[k] be the cumulative sum of frequencies up to 
%     the k-th sorted value. Then F[N] is the same as the total number of 
%     observations.
%
%     Define AvgFreq = F[N] / NumBins, which is the ideal or target average
%     frequency per bin. Similarly, the target cumulative frequency for the
%     first n bins is n*AvgFreq. The implementation of the equal-frequency
%     algorithm attempts to match the cumulative frequency up to the n-th
%     bin, and sets the n-th cut point by finding the index k such that the
%     distance abs(F[k] - n*AvgFreq) is minimized.
%
%     If a single value contains too many observations, equal-frequency
%     bins are not possible, and the above rule yields less than NumBins
%     total bins. In that case, the algorithm determines NumBins bins by
%     breaking up bins, in the order in which the bins were constructed.
%
%     The preprocessing of categorical predictors consists in sorting the 
%     categories according to the 'SortCategories' criterion (the default 
%     is to sort by odds in increasing order). Sorting is not applied to 
%     ordinal predictors. See Sort Categories in the Notes section below.
%
%   o Equal Width: Unsupervised algorithm that divides the range of values 
%     in the domain of the predictor variable into a predetermined number 
%     of bins of "equal width." For numeric data, the width is measured as 
%     the distance between bin edges. For categorical data, width is measured 
%     as the number of categories within a bin.For numeric predictors, if 
%     MinVal and MaxVal are the minimum and maximum data values in the 
%     sample, we let
%
%            Width = (MaxVal - MinVal)/NumBins
%
%     and the CutPoints are set to MinVal + Width, MinVal + 2*Width,
%     ... MaxVal - Width.
%
%     For categorical, if we have NumCats number of original categories, we
%     let
%
%            Width = NumCats / NumBins,
%
%     and set cut-point indices to one plus the rounded values of Width,
%     2*Width, ..., NumCats - Width.
%
%     The preprocessing of categorical predictors consists in sorting the 
%     categories according to the 'SortCategories' criterion (the default 
%     is to sort by odds in increasing order). Sorting is not applied to 
%     ordinal predictors. See Sort Categories in the Notes section below.
%
% Notes:
%
%  Sort Categories: As a preprocessing step for categorical data, 'Monotone', 
%  'EqualFrequency', and 'EqualWidth' support the 'SortCategories' input. 
%  This serves the purpose of reordering the categories before applying the 
%  main algorithm. The default sorting criterion is to sort by 'Odds'. For 
%  example, suppose that the data originally looks like this:
%
%           Bin	         Good	 Bad 	Odds
%           'Home Owner'  365	 177   2.062
%           'Tenant'	  307	 167   1.838
%           'Other'	      131	  53   2.472
%
%  After the preprocessing step, the rows would be sorted by 'Odds' and the
%  table looks like this:
% 
%           Bin	        Good	Bad     Odds
%           'Tenant'	 307	167	   1.838
%           'Home Owner' 365	177	   2.062
%           'Other'	     131	53	   2.472
%
%  The three algorithms only merge adjacent bins, so the initial order of 
%  the categories makes a difference for the final binning. The 'None' 
%  option for 'SortCategories' would leave the original table unchanged. 
%  For a description of the sorting criteria supported, see the description 
%  of the AlgorithmOptions option for 'SortCategories'.
% 
%  Upon the construction of a scorecard, the initial order of the categories, 
%  before any algorithm or any binning modifications are applied, is the 
%  order shown in the first output of BININFO. If the bins have been 
%  modified (either manually with modifybins or automatically with 
%  autobinning), use the optional output (cg,'category grouping') from 
%  BININFO to get the current order of the categories.
% 
%  The 'SortCategories' option has no effect on categorical predictors for 
%  which the 'Ordinal' parameter is set to true (see the 'Ordinal' input 
%  parameter in MATLAB categorical arrays for categorical. Ordinal data 
%  has a natural order, which is honored in the preprocessing step of the 
%  algorithms by leaving the order of the categories unchanged. Only 
%  categorical predictors whose 'Ordinal' parameter is false (default
%  option) are subject to reordering of categories according to the 
%  'SortCategories' criterion.
%
% References:
%
%   Refaat, M., "Data Preparation for Data Mining Using SAS," Morgan
%     Kaufmann, 2006.
%
%   Refaat, M., "Credit Risk Scorecards," lulu.com, 2011.
%
%   Anderson, R. "The Credit Scoring Toolkit," Oxford University Press,
%   2007.
%
%   Thomas, L., et al., "Credit Scoring and Its Applications", Society for 
%   Industrial and Applied Mathematics, 2002.
%
% See also BINDATA, BININFO, MODIFYBINS, PLOTBINS.
%

% Copyright 2014-2015 The MathWorks, Inc.

    SupportedAlgorithms = {'Monotone','EqualFrequency','EqualWidth'};

    parser = inputParser;

    parser.addOptional('PredictorVar',obj.PredictorVars,@(x)ischar(x)|iscell(x));
    parser.addParameter('Algorithm','Monotone',@(x)all(ismember(lower(x),...
                       lower(SupportedAlgorithms))));
    parser.addParameter('AlgorithmOptions',{},@(x)iscell(x));
    parser.addParameter('Display','off',@(x)ismember(lower(x),{'on','off'}));
    
    
    if ~isempty(varargin)
        % Locate where the string 'AlgorithmOptions' is in the input arguments
        % (when applicable)
        logicalIndex = cellfun(@(c)ischar(c)&strcmpi(c,'AlgorithmOptions'),...
            varargin,'UniformOutput',false);
        
        if iscell(logicalIndex)
            NumElems = cellfun(@numel,logicalIndex);
            if NumElems(1) > 1 % First input is a cell array of predictors
                logicalIndex = [logicalIndex{:}];
                logicalIndex(1:NumElems(1)-1) = [];
            else
                logicalIndex = [logicalIndex{:}];
            end
        end
        
        Ind = find(logicalIndex);

        if ~isempty(Ind)
            if numel(varargin) == Ind
                % Case where 'AlgorithmOptions' is the _last_ input string 
                % (entered by mistake)
                error(message('finance:creditscorecard:autobinning:IncompatibleNameValuePair'));

            else
                if ~iscell(varargin{Ind+1})
                    % Case where the parameter 'AlgorithmOptions' is passed,
                    % but its value is not a cell. We need to run this test 
                    % before we parse the input arguments to avoid mis-
                    % handling on the input arguments.
                    error(message('finance:creditscorecard:autobinning:InvalidAlgorithmOptions'));
                end
            end
        end

        if mod(numel(varargin),2) == 0
            % Even number of inputs means name-value pair parameters were 
            % entered. 
            % This is the case where no predictor name is passed. By default,  
            % the predictors in the object's property 'PredictorVars' are 
            % to be binned.
            % In case the first input is a predictor variable, an error is 
            % thrown.

            if any(ismember(varargin{1},obj.PredictorVars))
                error(message('finance:creditscorecard:autobinning:IncompatibleNameValuePair'));
            end

            inputArgs = [{obj.PredictorVars} varargin];

        else
            % Odd number of input arguments means that predictor variables 
            % were passed in, with potentially settings for the binning 
            % algorithm and algorithm options.
            if ~all(ismember(varargin{1},obj.PredictorVars))
                error(message('finance:creditscorecard:autobinning:InvalidPredictors'));
            end

            inputArgs = varargin;

        end
        
    else
        inputArgs = varargin;
        
    end
    
    parser.parse(inputArgs{:});

    PredictorVar     = parser.Results.PredictorVar;
    Algorithm        = parser.Results.Algorithm;
    AlgorithmOptions = parser.Results.AlgorithmOptions;
    Display          = parser.Results.Display;

    PredictorVar = obj.PreprocessPredictorVar(PredictorVar);

    %try
        [Algorithm,AlgorithmOptions] = obj.PreprocessAlgorithm(Algorithm,...
                                   AlgorithmOptions,PredictorVar);
    %catch MException
        %throw(MException);
    %end
    
    % If display is on, convert algorithm name to mixed case
    if strcmpi(Display,'on')
        [~,ind] = ismember(lower(Algorithm),lower(SupportedAlgorithms));
        Algorithm = SupportedAlgorithms(ind);
    end
    
    % Modifying the bins changes the predictor values of the logistic
    % model. Discard model information here, in case it was set before, to
    % force re-fitting the logistic model with the new binning.
    obj.ModelVars  = {};
    obj.ModelCoeff = [];
    
    % Format algorithm name into proper string for display
    AlgoName = Algorithm{1};
    [iStart,iEnd] = regexpi(AlgoName,'equal');
    if isempty(iStart)
        AlgoNameDisp = [upper(AlgoName(1)) lower(AlgoName(2:end))];
    else
        AlgoNameDisp = [upper(AlgoName(1)) lower(AlgoName(2:iEnd)) ' ' ...
            upper(AlgoName(iEnd+1)) lower(AlgoName(iEnd+2:end))];
    end
    
    for i = 1 :  numel(PredictorVar)
        if strcmpi(Display,'on')            
            fprintf('Binning ''%s'' with ''%s'' algorithm ... ',...
                PredictorVar{i},Algorithm{i});
        end
        %try
            %zhenxiao
            %obj = obj.BinDataByAlgorithm(PredictorVar(i),Algorithm,AlgorithmOptions);
            obj = obj.BinDataByAlgorithmByInd(i,Algorithm,AlgorithmOptions);
            obj.LatestBinning.(PredictorVar{i}) = ['Automatic / ' AlgoNameDisp];
        %catch MException
            %throw(MException);
        %end
        if strcmpi(Display,'on')
            fprintf('Done\n');
        end
    end

end
    
