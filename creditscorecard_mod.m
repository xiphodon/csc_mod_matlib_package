classdef (Sealed) creditscorecard_mod
% CREDITSCORECARD Create a creditscorecard object.
% 
% Create a CREDITSCORECARD object to perform manual and automaticn
% binning, fit a model using Weight Of Evidence data (WOE), determine
% scorecard points and score data.
%
% The typical workflow of using a creditscorecard class is as follows:
%
%   1) Create a creditscorecard object from a table. The toolbox includes
%   sample data from Refaat (2011).
%
%     Example:
%
%       load CreditCardData
%       sc = creditscorecard(data);
%       sc = creditscorecard(data,param1,val1,...);
%
%   2) Perform manual binning on a given predictor variable (see
%      MODIFYBINS).
%
%   3) Inquire about the bin information (see BININFO). Inquiring about
%      the bin information can be made anytime during the workflow. Bin
%      information includes Weight Of Evidence (WOE) and other measures.
%
%   4) The user can also plot the binned data (see PLOTBINS) and observe
%      how they change dynamically, as the binning is manually or
%      automatically modified (see MODIFYBINS, AUTOBINNING).
%
%   5) Alternatively, perform automatic binning using specific
%      algorithms (see AUTOBINNING).
%
%   6) Fit a logistic regression model to the binned data (see
%      FITMODEL). Alternatively, get WOE transformed data (see BINDATA),
%      fit a linear logistic regression model outside CREDITSCORECARD,
%      and set the model coefficients to CREDITSCORECARD (see SETMODEL).
%
%   7) Display and format creditscorecard points (see DISPLAYPOINTS,
%      FORMATPOINTS).
%
%   8) Score data (see SCORE).
%
%   9) Compute the probabilities of default (see PROBDEFAULT).
%
%  10) Validate the fitting model (see VALIDATEMODEL).
% 
%
% Input Arguments:
%     
%   data - Table containing numeric and/or categorical predictors.
%       Predictors with logical or string values are also accepted.
%          
% Optional Input Arguments and Name/Value Pairs:
%       
%   'GoodLabel'    Sets which of the two possible values in the response variable 
%                  correspond to 'Good' observations. By default this is the 
%                  value that has the highest number of observations in the
%                  response variable.
%
%   'IDVar'        String. Name of the variable used as ID for the
%                  observations. It is provided as a convenience to remove 
%                  the column corresponding to this variable from the predictor 
%                  variables. By default, it is empty.
%
%   'ResponseVar'  String. Name of the response variable. It must be binary. 
%                  By default, this is the last column of the data input.
%
%
% Creditscorecard public properties:
%
%   IDVar          Name of the variable used as ID or tag for the observations. 
%                  This property can also be set as an optional parameter 
%                  in the constructor.
%
%   PredictorVars  Cell array of predictor variables.
%
%
% Creditscorecard public get access properties:
%
%   GoodLabel      Sets which of the two possible values in the response
%                  variable correspond to 'Good' observations. This
%                  property can only be set with an optional parameter in 
%                  the constructor.
%
%   VarNames       Cell array of predictor variable names. These come
%                  directly from the data input in the constructor.
%
%   ResponseVar    Name of the response variable. This property can only be
%                  set with an optional parameter in the constructor.
%
%   Data           Table object, containing the data used to create the
%                  creditscorecard.
%
%
% Creditscorecard methods:
%
%   bininfo        Returns statistics about the data as well as the cutpoints 
%                  or category groupings.
%
%   modifybins     Performs manual binning on the data, by changing the cut 
%                  points or category groupings.
%
%   autobinning    Performs automatic binning by applying one of the supported 
%                  algorithms.
%
%   plotbins       Plots a histogram of the binned data, as well as the WOE 
%                  curve, with information on the bin count.
%
%   bindata        Returns a table of binned predictors. This method is
%                  used within FITMODEL to bin the WOE data. 
%
%   fitmodel       Fits a linear regression model to the WOE data.
%
%   setmodel       An alternative to FITMODEL. Directly sets the model predictor 
%                  names and coefficients.
%
%   formatpoints   Formats the creditscorecard points.
%
%   displaypoints  Returns the creditscorecard points information per 
%                  predictor, per bin.
%
%   score          Determines scores (and points) for given data.
%
%   probdefault    Compute the probabilities of default for given data.
%
%   validatemodel  Computes statististics (CAP, ROC, KS) to validate the
%                  fitting model.
%
%
% References:
%
%   Anderson, R., "The Credit Scoring Toolkit," Oxford University
%      Press, 2007.
%
%   Refaat, M., "Credit Risk Scorecards," lulu.com, 2011.
%
%   See also MODIFYBINS, BININFO, PLOTBINS, AUTOBINNING, FITMODEL, 
%   DISPLAYPOINTS, SCORE, PROBDEFAULT, VALIDATEMODEL.

%   Copyright 2014-2016 The MathWorks, Inc.
    
    
properties (Access = public)
    % Data properties
    RawData
    ResponseVarMap
    ResponseOrder
    
    % Binning properties
    BinContainers = {}; 
    LatestBinning = struct();

    % For model fitting
    ModelVars  = {};
    ModelCoeff = [];
    
    % Plot properties
    BinText   = {};
    HasWOE    = false;
    
    % To modify and scale points and score
    PointsFormat = struct('BasePoints'    , false,...
                          'Round'         , 'None',...
                          'Missing'       , 'NoScore',...
                          'ShiftAndSlope' , [0 1],...
                          'WorstAndBestScores', [],...
                          'PointsOddsAndPDO', []);
    
end

properties(Access = private)
    PrivateIDVar = '';
    PrivatePredictorVars = {};
end

properties(SetAccess = private, GetAccess = public)
    GoodLabel          % The GoodLabel value will be mapped to 1. 
    ResponseVar = '';  % Name of the response variable.
    VarNames           % Cell array of predictor variable names.
    NumericPredictors = {};
    CategoricalPredictors = {};
end

properties(Dependent)
    IDVar = '';
    PredictorVars
end

properties(Dependent, SetAccess = private, GetAccess = public)
    Data
end

methods
    % Creditscorecard Constructor
    
    function obj = creditscorecard_mod(data,varargin)
        % CREDITSCORECARD Create a creditscorecard object.
        %
        % Syntax:
        %
        %   sc = creditscorecard(data)
        %   sc = creditscorecard(data,param1,val1,...)
        %
        % Description:
        %
        %   Create a creditscorecard object by specifying training data.
        %
        % Input Arguments:
        %
        %   data - A table containing training data. Each column must be 
        %       either a numeric or a categorical array. It must contain 
        %       a binary (i.e. two-valued) response variable.
        %
        %
        % Optional Input Parameter Name/Value Pairs:
        %
        %   'ResponseVar'   A string containing the name of the column that
        %                   corresponds to the response variable. If not specified,
        %                   it is set to the last column of the data input.
        %
        %   'GoodLabel'     A numeric scalar (if response variable is numeric) 
        %                   or a  string (if response is categorical) indicating 
        %                   which of the two possible response values corresponds 
        %                   to 'Good'. If not specified, it is set to the 
        %                   response value with the highest count.
        %
        %   'IDVar'         A string containing the name of the variable
        %                   which is treated as an ID. 
        % 
        %   'PredictorVars' A cell of strings containing the names of the
        %                   variables to be considered later for binning
        %                   the data and fitting a model.
        %
        % Output:
        %
        %   sc - A creditscorecard object.
        %
        
        if nargin == 0
            error(message('finance:creditscorecard:creditscorecard:MissingInputArgument'));
        end
        
        if ~istable(data)
            error(message('finance:creditscorecard:creditscorecard:DataFormat'));
        end
        
        obj.VarNames = data.Properties.VariableNames;
                               
        parser = inputParser;

        parser.addParameter('ResponseVar','',@(x)ischar(x)|iscell(x));
        parser.addParameter('GoodLabel','',@(x)(isnumeric(x)||ischar(x)||islogical(x)));
        parser.addParameter('IDVar','',@(x)(ischar(x)||iscell(x)));
        parser.addParameter('PredictorVars',{},@(x)(ischar(x)||iscell(x)));
        
        parser.parse(varargin{:});
        
        obj.ResponseVar = obj.PreprocessResponseVar(parser.Results.ResponseVar,...
                                                    obj.VarNames);
                                                
        % Validate predictor names and IDVar
        if ~isempty(parser.Results.IDVar) && ~ismember(parser.Results.IDVar,obj.VarNames)
             error(message('finance:creditscorecard:creditscorecard:InvalidIDVar'));
        end
        
        if any(~ismember(parser.Results.PredictorVars,obj.VarNames))
            error(message('finance:creditscorecard:creditscorecard:InvalidPredictorName'));
        end

        
        % Validate data type for predictors
        obj.RawData = obj.PreprocessPredictorData(data);
        
        if isempty(parser.Results.PredictorVars)
            obj.PredictorVars = obj.VarNames;
        else
            obj.PredictorVars = parser.Results.PredictorVars;
        end
        
        obj.IDVar = parser.Results.IDVar;
        
        obj.GoodLabel = parser.Results.GoodLabel;
        
        ResponseValues = obj.validateResponseValues(obj.RawData.(obj.ResponseVar));
        
        if ~isempty(obj.GoodLabel)
            obj = obj.validateGoodLabel(ResponseValues);
        end
        
        [obj,mappedRespVar] = obj.MapResponseVar(obj.RawData.(obj.ResponseVar),...
                                ResponseValues);

        if ~isnumeric(obj.GoodLabel)
            obj.ResponseOrder = {obj.GoodLabel char(setdiff(ResponseValues,...
                                 obj.GoodLabel))};
        else
            obj.ResponseOrder = [obj.GoodLabel setdiff(ResponseValues,...
                                 obj.GoodLabel)];
        end
        
        obj.ResponseVarMap = mappedRespVar;
        
        % Initialize the LatestBinning property. Upon creation, the property
        % is set to "Original Data" for all predictors
        
        for i = 1 : length(obj.PredictorVars)
            Predictor = obj.PredictorVars{i};
            obj.LatestBinning.(Predictor) = 'Original Data';
        end
        
    end

    function value = get.VarNames(obj)
        value = obj.VarNames;
    end    
    
    function obj = set.IDVar(obj,idvar)
        if ~isempty(idvar) 
            if ~ischar(idvar)
                if iscell(idvar)
                    if numel(idvar) > 1
                        error(message('finance:creditscorecard:creditscorecard:InvalidIDVarNum'));
                        
                    else
                        if ~ismember(idvar{:},obj.VarNames) || ...
                                        strcmpi(idvar{:},obj.ResponseVar)
                            error(message('finance:creditscorecard:creditscorecard:InvalidIDVar'));
                            
                        else
                            obj.PrivateIDVar = idvar{:};
                        end
                    end
                    
                else
                    error(message('finance:creditscorecard:creditscorecard:InvalidIDVarType'));
                    
                end
                
            else
                
                if ~ismember(idvar,obj.VarNames) || strcmpi(idvar,obj.ResponseVar)
                    error(message('finance:creditscorecard:creditscorecard:InvalidIDVar'));

                else
                    obj.PrivateIDVar = idvar;
                end
                
            end
        else
            obj.PrivateIDVar = '';
        end
        
    end
    
    function idvar = get.IDVar(obj)
        idvar = obj.PrivateIDVar;
    end
    
    function obj = set.PredictorVars(obj,predictorVars)
        if any(~ismember(predictorVars,obj.VarNames))
            error(message('finance:creditscorecard:creditscorecard:InvalidPredictorName'));
        end
        obj.PrivatePredictorVars = predictorVars;
    end
    
    function predictorVars = get.PredictorVars(obj)
        if isempty(obj.PrivatePredictorVars)
            obj.PrivatePredictorVars = obj.VarNames;
        end
        predictorVars = setdiff(obj.PrivatePredictorVars,{obj.IDVar,obj.ResponseVar},'stable');
    end  

    function NumPredictors = get.NumericPredictors(obj)
        NumPredictors = cell(size(obj.PredictorVars));
        for i = 1 : length(NumPredictors)
            data = obj.RawData.(obj.PredictorVars{i});
            if isnumeric(data) 
                NumPredictors{i} = obj.PredictorVars{i};
            end
        end
        NumPredictors(cellfun('isempty',NumPredictors)) = [];
    end  

    function CatPredictors = get.CategoricalPredictors(obj)
        CatPredictors = cell(size(obj.PredictorVars));
        for i = 1 : length(CatPredictors)
            data = obj.RawData.(obj.PredictorVars{i});
            if iscategorical(data)
                CatPredictors{i} = obj.PredictorVars{i};
            end
        end
        CatPredictors(cellfun('isempty',CatPredictors)) = [];
    end  
    
    function data = get.Data(obj)
        data = obj.RawData;
    end
    
end

%zhenxiao
methods (Access = public, Hidden = false)
    %% Process predictors and response variable
    
    function ResponseVar = PreprocessResponseVar(~,ResponseVar,VarNames)
        % Preprocess response variable. If the user does not explicitly
        % specify which variable will be used as the response variable,then
        % by default, the last column of the table is used.
        
        NumVars = numel(VarNames);
        if isempty(ResponseVar)
           % The response variable is the last one by default
           RespInd = NumVars;
        else
           if isnumeric(ResponseVar)
              if any(~ismember(ResponseVar,1:NumVars))||numel(ResponseVar)~=1
                 error(message('finance:creditscorecard:creditscorecard:InvalidResponseIndex')); 
              end
              RespInd = ResponseVar;
           else
              [ok,Ind] = ismember(ResponseVar,VarNames);
              if ~ok
                 error(message('finance:creditscorecard:creditscorecard:InvalidResponseName')); 
              end
              RespInd = Ind;
           end
        end
        ResponseVar = VarNames{RespInd};
    end
    
    
    function PredictorVar = PreprocessPredictorVar(obj,PredictorVar)                      
        % Preprocess the predictor variable names. Check whether or not they 
        % belong to the VarNames property of the creditscorecard object. 
        % Otherwise, issue an error.
        
        NumVars = numel(obj.VarNames);
        [~,RespInd] = ismember(obj.ResponseVar,obj.VarNames);
        [~,IDInd]   = ismember(obj.IDVar,obj.VarNames);
        VarsInd     = setdiff(1:NumVars,[IDInd RespInd]);
        VarNamesFiltered = obj.VarNames(VarsInd);
        
        if isempty(PredictorVar)
            PredictorInd = VarsInd;   
            PredictorVar = obj.VarNames(PredictorInd); 
        else
            if ismember(PredictorVar,{obj.IDVar obj.ResponseVar})
                error(message('finance:creditscorecard:creditscorecard:InvalidPredictorName'))
            end
            
            [ok,Ind] = ismember(PredictorVar,VarNamesFiltered);
            if any(~ok)
                error(message('finance:creditscorecard:creditscorecard:InvalidPredictorName')); 
            end
            
            PredictorInd = Ind;
            PredictorVar = VarNamesFiltered(PredictorInd); 
        end
        
        if isempty(PredictorVar)
            error(message('finance:creditscorecard:creditscorecard:EmptyPredictor'));
        end
        
    end
    
    
    function data = PreprocessPredictorData(~,data)
        % Helper method that validates the data type for predictors. Accepted
        % types are numeric and categorical. If logical, cellstr and character
        % arrays are passed, they are converted to numeric for the former, 
        % and categorical for the latter two. 
        %
        % Input parameters:
        %
        %   data - Table of predictors
        %
        % Outputs:
        %
        %   data - Updated table to reflect the conversion explained above.
        %
        
        VarsNames = data.Properties.VariableNames;
        NumVars   = numel(VarsNames);
        
        for i = 1:NumVars
            if ~iscategorical(data.(VarsNames{i})) && ~isnumeric(data.(VarsNames{i}))
                
                if islogical(data.(VarsNames{i}))
                    % Convert to double
                    data.(VarsNames{i}) = double(data.(VarsNames{i}));
                    
                elseif ischar(data.(VarsNames{i}))
                    % Convert to categorical
                    data.(VarsNames{i}) = categorical(cellstr(data.(VarsNames{i})));
                    
                elseif iscellstr(data.(VarsNames{i}))
                    % Convert to categorical
                    data.(VarsNames{i}) = categorical(data.(VarsNames{i}));
                    
                elseif isstring(data.(VarsNames{i}))
                    % Convert to categorical
                    data.(VarsNames{i}) = categorical(data.(VarsNames{i}));
       
                else
                    error(message('finance:creditscorecard:creditscorecard:DataType'));
                    
                end
            end
        end
        
    end
    
    
    function ResponseValues = validateResponseValues(~,ResponseData)
        
        ResponseValues = unique(ResponseData);

        ResponseValues(isnan(double(ResponseValues))) = [];
        
        if numel(ResponseValues) ~= 2
            error(message('finance:creditscorecard:creditscorecard:ResponseValues'));   
        end
        
    end
    
    
    function obj = validateGoodLabel(obj,ResponseValues)
        
        if islogical(obj.GoodLabel)
            obj.GoodLabel = double(obj.GoodLabel);
            
        elseif ischar(obj.GoodLabel)
            obj.GoodLabel = categorical(cellstr(obj.GoodLabel));
            
        end
        
        if ~strcmpi(class(obj.GoodLabel),class(ResponseValues))
            error(message('finance:creditscorecard:creditscorecard:InvalidGoodLabel'));
        end
        
        if ~ismember(obj.GoodLabel,ResponseValues)
            error(message('finance:creditscorecard:creditscorecard:InvalidGoodLabel'));
        end
        
    end
    
    
    function [obj,mappedResponseVar] = MapResponseVar(obj,ResponseData,...
                                                ResponseValues)
        % ResponseData: Array containing the reponse values for given
        % observations.
        
        RespVar = ResponseData; 
        RespNaNInd = isnan(double(ResponseData));
        
        if ~isempty(obj.GoodLabel)       
            if isnumeric(obj.GoodLabel)               
               RespVar = double(RespVar == obj.GoodLabel);
               
            else
                GoodInd = (RespVar == obj.GoodLabel);
                RespVar = double(RespVar)-1;
                RespVar(GoodInd)  = 1;
                RespVar(~GoodInd) = 0;
                
                % Convert back the GoodLabel, for display purposes
                obj.GoodLabel = char(obj.GoodLabel);
            end

        else
            % No 'GoodLabel' is supplied. Do a count of the predominant
            % response variable value and map the corresponding entries 
            % to '1'.
            if isnumeric(ResponseValues) 
                s1 = sum(ResponseData == ResponseValues(1));
                s2 = sum(ResponseData == ResponseValues(2));
                if s1 > s2
                    RespVar = double(ResponseData == ResponseValues(1));
                    obj.GoodLabel = ResponseValues(1);
                    
                else
                    RespVar = double(ResponseData == ResponseValues(2));
                    obj.GoodLabel = ResponseValues(2);
                    
                end  
            else
                
                RespVar = double(ResponseData)-1;
                s1 = sum(RespVar == 0); % ResponseValues(1) = category1 is mapped to 0
                s2 = sum(RespVar == 1); % ResponseValues(2) = category2 is mapped to 1
                if s1 > s2
                    RespVar = double(ResponseData == ResponseValues(1));
                    obj.GoodLabel = ResponseValues(1);
                    
                else
                    RespVar = double(ResponseData == ResponseValues(2));
                    obj.GoodLabel = ResponseValues(2);
                end   
                
                % Cast to character for display purposes
                obj.GoodLabel = char(obj.GoodLabel);
            end
        end
        
        RespVar(RespNaNInd) = NaN;
        mappedResponseVar = RespVar;
    end
    
    function ValidateDataType(obj,dbc,DataStruct)      
        % Called by MODIFYBINS. Validate that the data type are consistent 
        % with the groupings. If the predictor variable is numeric, then
        % the user must provide a numerical array of cut points. If the
        % predictor variable is categorical, then the user must provide a
        % table with category names and corresponding bin numbers.
              
        if ismember(DataStruct.PredictorVar,{obj.ResponseVar obj.IDVar})
            error(message('finance:creditscorecard:creditscorecard:InvalidPredictorToBin'));
        end
        
        if isnumeric(dbc.PredictorData)
            if ~isempty(DataStruct.CutPoints) && ~isnumeric(DataStruct.CutPoints)
                error(message('finance:creditscorecard:creditscorecard:DataTypeMismatch')); 
            end
            
            if ~isempty(DataStruct.CatGrouping)
                error(message('finance:creditscorecard:creditscorecard:DataTypeMismatch')); 
            end
        end
        
        if iscategorical(dbc.PredictorData)
            if ~isempty(DataStruct.CatGrouping) && ~istable(DataStruct.CatGrouping)
                error(message('finance:creditscorecard:creditscorecard:DataTypeMismatch')); 
            end
            
            if ~isempty(DataStruct.CutPoints)
                error(message('finance:creditscorecard:creditscorecard:DataTypeMismatch')); 
            end
            
            if ~isempty(DataStruct.MinValue) || ~isempty(DataStruct.MaxValue)
                error(message('finance:creditscorecard:creditscorecard:DataTypeMismatchMinMax')); 
            end
            
        end
    end
    
    
    %% Process algorithms and algorithm options
    
    function [Algo,AlgoOptions] = PreprocessAlgorithm(obj,Algo,AlgoOptions,...
                                  PredictorVar)                      
        % Pre-process algorithm name and options. Uses default
        % 'EqualFrequency' algorithm with five bins.
        
        NumPredictors     = numel(PredictorVar);
        SingleAlgo        = false;
        SingleAlgoOptions = false;
        
        % Handling of algorithms
        if ischar(Algo)
            Algo       = repmat({Algo},NumPredictors,1);
            SingleAlgo = true;
        end
        
        % Handling of algorithm options
        if isempty(AlgoOptions)
            AlgoOptions = obj.SetDefaultAlgorithmOptions(Algo);
        end
        
        % First, scalar expansion
        if ~iscell(AlgoOptions{1}) ||...
                (iscell(AlgoOptions{1}) && numel(AlgoOptions) == 1)
            AlgoOptions = repmat({AlgoOptions},NumPredictors,1);
            SingleAlgoOptions = true;
        end
        
        % Check size of inputs
        if numel(AlgoOptions) ~= NumPredictors
            error(message('finance:creditscorecard:creditscorecard:AlgorithmOptionsMismatch')); 
        end
        
        % Construct algorithms and set options to catch errors in the input 
        % before starting the binning process
        if SingleAlgo && SingleAlgoOptions
            %try
                [~] = obj.CreateBinningAlgorithm(Algo{1},AlgoOptions{1});
            %catch
                %error(message('finance:creditscorecard:creditscorecard:InvalidAlgorithmSpecs')); 
            %end
        else
            for i = 1 : numel(Algo)
                try
                    [~] = obj.CreateBinningAlgorithm(Algo{i},AlgoOptions{i});
                catch
                    error(message('finance:creditscorecard:creditscorecard:InvalidAlgorithmSpecs')); 
                end
            end
        end
        
    end
 
    
    function obj = BinDataByAlgorithm(obj,PredictorVars,Algorithm,AlgorithmOptions)
        % Apply the algorithm to bin the given predictor variables. If
        % more than one predictor variable is passed, the same algorithm and
        % algorithm options are used for all.
         
        NumPredictors = numel(PredictorVars);
        
        for i = 1 : NumPredictors
            % Create or retrieve data container
            [obj,Ind] = obj.getContainerIndex(PredictorVars{i});
            
            % Create algorithm and set options
            algo = obj.CreateBinningAlgorithm(Algorithm{i},...
                   AlgorithmOptions{i});
            
            % Run algorithm
            obj.BinContainers{Ind} = algo.runAlgorithm(obj.BinContainers{Ind});
            
        end
        
    end
    
    %zhenxiao
    function obj = BinDataByAlgorithmByInd(obj,Ind,Algorithm,AlgorithmOptions)
       obj = obj.getContainerByInd(Ind);
            
       % Create algorithm and set options
       algo = obj.CreateBinningAlgorithm(Algorithm{Ind},...
       AlgorithmOptions{Ind});
            
       % Run algorithm
        obj.BinContainers{Ind} = algo.runAlgorithm(obj.BinContainers{Ind});
            
    end
        
   
    function algooptions = SetDefaultAlgorithmOptions(obj,AlgoType)        
        % If no algorithm options are passed, default values and parameters
        % are used.
        
        AlgoMap     = obj.GetBinningAlgorithmMap; 
        [~,ind]     = ismember(lower(AlgoType),AlgoMap(:,1));
        algooptions = AlgoMap{ind,3};
        
    end
    
    
    function algo = CreateBinningAlgorithm(obj,AlgoType,AlgoOptions)
        
        AlgoMap = obj.GetBinningAlgorithmMap(); 
        [~,ind] = ismember(lower(AlgoType),AlgoMap(:,1)); 
        algo    = AlgoMap{ind,2}();  
        algo    = algo.setOptions(AlgoOptions{:});
        
    end  
    
    
    function AlgoMap = GetBinningAlgorithmMap(~)
        % The following are the default options for each algorithm
        
        AlgoMap = {...
            'equalfrequency',@EqualFrequency_mod,{'NumBins',5};
            'equalwidth',@EqualWidth_mod,{'NumBins',5};
            'monotone',@Monotone_mod,{'InitialNumBins',10}};
    end
    
    
    %% Binning container helper methods
    
    function dbc = CreateBinningContainer(~,PredictorData,ResponseData,varargin)
        
        if iscategorical(PredictorData)
            dbc = Categorical_mod(PredictorData,...
                ResponseData,varargin{:});
        else
            dbc = Numeric_mod(PredictorData,...
                ResponseData,varargin{:});
        end
    end
    
    
    function [obj,Ind] = getContainerIndex(obj,varargin)
        % Checks if a container for the given predictor variable already 
        % exists in the cell array BinContainers, property of the creditscorecard 
        % object. If the container does not exist, then it is created. This
        % helper method is called by MODIFYBINS and AUTOBINNING.
        %
        % varargin{1} - string containg the name of the predictor to bin
        % varargin{2} - If any, over-writes the binning container with the
        %               one passed in as argument.
        
        PredictorName = varargin{1};
        
        if iscell(PredictorName)
            PredictorName = PredictorName{:};
        end
        
        if isempty(obj.BinContainers)
            if nargin > 2
                dbc = varargin{2};
            else
                dbc = obj.CreateBinningContainer(obj.RawData.(PredictorName),...
                    obj.RawData.(obj.ResponseVar),'PredictorVar',PredictorName,...
                    'ResponseVar',obj.ResponseVar,'ResponseOrder',...
                    obj.ResponseOrder);
                
            end
            obj.BinContainers = {dbc};
            Ind = 1;
            
        else
            [ok,Ind] = ismember(PredictorName,cellfun(@(c)c.PredictorVar,...
                obj.BinContainers,'UniformOutput',false));
            if ok
                % Only need to update the 'BinContainers' stack if the dbc
                % is passed. Otherwise, we only need 'Ind' above .
                if nargin > 2
                    obj.BinContainers{Ind} = varargin{2};
                end
                
            else
                % The dbc is not in the 'BinContainers' stack
                if nargin > 2
                    dbc = varargin{2};
                else
                    % This is the case where we are binning a new predictor
                    % variable, not already stored in the 'BinContainer'
                    % property of the creditscorecard.
                    dbc = obj.CreateBinningContainer(obj.RawData.(PredictorName),...
                        obj.RawData.(obj.ResponseVar),'PredictorVar',PredictorName,...
                        'ResponseVar',obj.ResponseVar,'ResponseOrder',...
                        obj.ResponseOrder);
                    
                end
                obj.BinContainers{end+1,:} = dbc;
                Ind = numel(obj.BinContainers);
                
            end
            
        end
    end
    
    
    %zhenxiao
    function obj = getContainerByInd(obj,Ind)
        PredictorName = obj.PredictorVars{Ind};
        if isempty(obj.BinContainers)
            obj.BinContainers = cell(size(obj.PredictorVars, 2), 1);
        end
        if isempty(obj.BinContainers{Ind, :})
            obj.BinContainers{Ind, :} = obj.CreateBinningContainer(obj.RawData.(PredictorName),...
                    obj.RawData.(obj.ResponseVar),'PredictorVar',PredictorName,...
                    'ResponseVar',obj.ResponseVar,'ResponseOrder',...
                    obj.ResponseOrder);
        end
         
    end
    
    
    function obj = removeBinContainer(obj,PredictorName)
        % Remove bin containers to keep consistency in data types.
        
        [obj,Ind] = obj.getContainerIndex(PredictorName);
        obj.BinContainers(Ind) = [];
        
    end
  
    
    %% Graphics helper methods
    
    function hFig = formatBarPlot(~,FreqTable,dbc,Normalized,iswoe,haslegend)
        % Helper method called by PLOTBINS. Displays the bar plot and sets 
        % the color, labels and bincounts of the bins.
        
        hFig = figure;
        hAx  = axes('Parent',hFig);
        
        if dbc.HasResponseData
            if size(FreqTable,1) == 1
                hBar = bar([FreqTable;nan(size(FreqTable))],'Stacked',...
                    'Parent',hAx);
            else
                hBar = bar(FreqTable,'Stacked','Parent',hAx);
            end
            set(hBar(1),'FaceColor',[102 153 255]/255,'Tag',dbc.PredictorVar);
            set(hBar(2),'FaceColor',[255 102 102]/255,'Tag',dbc.PredictorVar);
            
        else
            hBar = bar(sum(FreqTable,2),'Parent',hAx);
            set(hBar,'FaceColor',[102 153 255]/255,'Tag',dbc.PredictorVar)
            
        end
       
        XTagLoc = get(hAx,'XTick');
        YTagLoc = sum(FreqTable,2);
        % In case unbinned data are used, the size of XTick might be
        % unreliable. Need to check size of XTagLoc.
        if numel(XTagLoc) ~= numel(YTagLoc)
            if any(isnan(get(hBar(1),'YData')))
                XTagLoc = 1;
            else
                XTagLoc = get(hBar(1),'XData')';
            end
            set(hAx,'XTick',XTagLoc);
        end
        
        yLims   = get(hAx,'YLim');
        set(hAx,'XTickLabel',dbc.getLabels);
        set(hAx,'Ylim',[yLims(1) 1.15*yLims(2)]);
        set(hAx,'Xlim',[0.5 dbc.NumBins+0.5]); 
        
        % ========= Set the Y-axis label =========
        set(get(hAx,'YLabel'),'String','Bin count');
        
        % ========= Format the bin count labels =========
        FreqLabels = internal.finance.binning.utils.FormatBinLabels(FreqTable,dbc.HasResponseData,dbc.NumBins,Normalized);
        
        if ~isempty(FreqLabels)
            if dbc.HasResponseData
                text(XTagLoc,YTagLoc,strcat(FreqLabels(:,1),FreqLabels(:,2),...
                    FreqLabels(:,3)),'FontWeight','Bold','HorizontalAlignment',...
                    'Center','VerticalAlignment','Bottom','Parent',hAx);
            else
                text(XTagLoc,YTagLoc,FreqLabels,'FontWeight','Bold',...
                    'HorizontalAlignment','Center','VerticalAlignment',...
                    'Bottom','Parent',hAx);
            end
        end
        
        % ========= Plot the WOE curve =========
        if strcmpi(iswoe,'on')
            WOE = internal.finance.binning.utils.WeightOfEvidenceByFreq(FreqTable);
            hold(hAx,'on')
            [AX,h1,h2] = plotyy(NaN,NaN,1:(dbc.NumBins),WOE,'plot','plot',...
                'Parent',hAx);
            set(AX(2),'Xlim',get(AX(1),'Xlim'));
            set(AX(2),'XTick',[],'YColor',[0 0 0],'Tag','WOE');
            set(get(AX(2),'YLabel'),'String','WOE');
            set(AX(1),'YColor',[0 0 0]);
            set(h2,'Color',[0 0 0],'Marker','o','MarkerFaceColor',...
                [0 0 0],'Linewidth',1);
            delete(h1)
            hold(hAx,'off');
            
        end
        
        % ========= Set the legend and title =========
        if dbc.HasResponseData && strcmpi(haslegend,'On')
            hl = legend('Good','Bad','location','SouthOutside','Orientation',...
                'Horizontal');
            set(hl,'Parent',hFig);
            updateAxesLegendPosition(hAx,hl);
        end
        
        title(hAx,dbc.PredictorVar,'Interpreter','None')
        
        function updateAxesLegendPosition(ha,hLegend)
            posA = get(ha,'Position');
            posL = get(hLegend,'Position');
            set(ha,'Position',[posA(1) 0.19 posA(3) 0.72]);
            set(hLegend,'Position',[posL(1) 0.05 posL(3) posL(4)]);
        end
        
    end
    
  
    %% Helper method to compute points
    
    function [S,Smin,Smax] = getUnscaledPoints(obj)
        % Description:
        %
        %  Helper method that computes unscaled points for the given
        %  creditsorecard object. Used by FORMATPOINTS to get the minimum 
        %  and maximum unscaled scores, and by GETFORMATTEDPOINTS to apply 
        %  scaling and pass scaled points to SCORE and DISPLAYPOINTS.
        %
        %  The computation of unscaled points is isolated because there are
        %  use cases that require unscaled points and scores directly. For
        %  example, FORMATPOINTS requires the unscaled minimum and maximum
        %  scores to determine scaling parameters.
        %
        % Output arguments:
        %
        %   S     Structure array whose field names are the predictor 
        %         variables used in the fitting model. The values of 
        %         'S' are tables with bin labels and unscaled points
        %         (without base points) corresponding to each bin. There
        %         is a separate field to report the unscaled base points,
        %         which, in the unscaled case, is simply the 'beta0' from
        %         the model coefficients, if available (zero otherwise).
        %
        %   Smin  Minimum possible unscaled score.
        %
        %   Smax  Maximum possible unscaled score.
        % 
        
        if isempty(obj.ModelVars) || isempty(obj.ModelCoeff)
            error(message('finance:creditscorecard:creditscorecard:EmptyModel'));
        end
        
        if (numel(obj.ModelCoeff) == 1 && strcmpi(obj.ModelVars{:},'(intercept)'))
            error(message('finance:creditscorecard:creditscorecard:InterceptOnlyModel'));
        end
        
        if ismember('(intercept)',lower(obj.ModelVars))
            % There are N coefficients, N-1 predictors and an intercept
            S.BasePoints = obj.ModelCoeff(1);
            modelVars   = obj.ModelVars(2:end);
            modelCoeffs = obj.ModelCoeff(2:end);
        else
            S.BasePoints = 0;
            modelVars   = obj.ModelVars;
            modelCoeffs = obj.ModelCoeff;
        end
        
        if ~isempty(modelVars)
           
            Smin = 0;
            Smax = 0;
        
            for i = 1 : numel(modelVars)

                [obj,Ind] = obj.getContainerIndex(modelVars{i});   
                Stats     = obj.BinContainers{Ind}.getStatistics('StatsList','woe');
                woedata   = Stats.woe;

                BinLabels = obj.BinContainers{Ind}.getLabels();

                % Exclude base points here, they are reported separately
                points    = modelCoeffs(i) * woedata(:);

                Smin = Smin + min(points);
                Smax = Smax + max(points);

                S.(modelVars{i}) = table(BinLabels(:),points,...
                    'VariableNames',{modelVars{i} 'Points'});

            end
            
            % Add unscaled base points to unscaled min and max scores
            Smin = S.BasePoints + Smin;
            Smax = S.BasePoints + Smax;
            
        else
            % The ModelCoeff property of the creditscorecard object 
            % contains only the 'intercept'. All model predictors  were 
            % dropped during the fitting. Throw an error.
            error(message('finance:creditscorecard:creditscorecard:EmptyModel'));
        end
        
    end
    
    
    function [S,Smin,Smax] = getFormattedPoints(obj)
       % Description:
       %
       %  This is where the scaling, and some formatting at the point
       %  level, takes place. This method gets unscaled points, and then applies
       %  scaling separation of base points (if requested), and rounding at
       %  the points level (if requested).
       %
       % Output arguments:
       %
       %   S     Structure array whose fieldnames are the predictor 
       %         variables used in the fitting model. The values of 
       %         'S' are tables with bin labels and scaled and rounded
       %         points (if requested).
       %         There is a separate field to report the scaled and
       %         rounded base points, if requested.
       %         If base points are not required to be reported
       %         separately, the base points are divided across the
       %         predictors points, and the BasePoints field is zero.
       %
       %   Smin  Minimum possible scaled score. This is naturally rounded
       %         if all points are rounded. Minimum is the lowest possible
       %         total score in the mathematical sense, independently of
       %         whether a low score means high risk or low risk.
       %
       %   Smax  Maximum possible scaled score. This is naturally rounded
       %         if all points are rounded. Maximum is the highest possible
       %         total score in the mathematical sense, independently of
       %         whether a low score means high risk or low risk.
       %

       [S,SUmin,SUmax] = obj.getUnscaledPoints;

       FieldNames = fields(S);
       PredNames = setdiff(FieldNames,'BasePoints','stable');
       NumPreds = length(PredNames);
       % Add field for predictor points when 'Missing' is set to 'ZeroWOE'
       S.ZeroWOEValue = 0;
       
       ShiftAndSlope = obj.PointsFormat.ShiftAndSlope;
       if isempty(ShiftAndSlope)
          if ~isempty(obj.PointsFormat.WorstAndBestScores)
             % Solve:
             %   WorstScore = Shift + Slope*SUmin
             %   BestScore  = Shift + Slope*SUmax
             ShiftAndSlope = [1 SUmin; 1 SUmax] \ obj.PointsFormat.WorstAndBestScores(:);
          elseif ~isempty(obj.PointsFormat.PointsOddsAndPDO)
             % Solve:
             %   P = Shift + Slope*log(Odds)
             %   P+PDO = Shift + Slope*log(2*Odds)
             P = obj.PointsFormat.PointsOddsAndPDO(1);
             Odds = obj.PointsFormat.PointsOddsAndPDO(2);
             PDO = obj.PointsFormat.PointsOddsAndPDO(3);
             ShiftAndSlope = [1 log(Odds); 1 log(2*Odds)] \ [P; P + PDO];
          end
       end
       
       A = ShiftAndSlope(1);
       B = ShiftAndSlope(2);
       
       % Compute points, depending on base points
       if obj.PointsFormat.BasePoints
          
          S.BasePoints = A + B*S.BasePoints;
          
          for i = 1 : NumPreds
             S.(PredNames{i}).Points = B* S.(PredNames{i}).Points;
          end
          
       else
          
          BasePointsPred = (A + B*S.BasePoints) / NumPreds;
          S.BasePoints = 0;
          S.ZeroWOEValue = BasePointsPred;
          for i = 1 : NumPreds
             S.(PredNames{i}).Points = BasePointsPred + B*S.(PredNames{i}).Points;
          end
          
       end

       % Round if all points must be rounded
       if strcmpi(obj.PointsFormat.Round,'allpoints')
          % Round all points, and reset minimum and maximum scores (rounding 
          % points may change total minimum and maximum).
          S.BasePoints = round(S.BasePoints);
          Smin = S.BasePoints;
          Smax = S.BasePoints;
          
          for i = 1 : NumPreds
             S.(PredNames{i}).Points = round(S.(PredNames{i}).Points);
             Smin = Smin + min(S.(PredNames{i}).Points);
             Smax = Smax + max(S.(PredNames{i}).Points);
          end
          
       end
       
       % Get scaled minimum and maximum scores
       % Note these are minimum and maximum values AFTER the scaling. Minimum
       % Score is always the minimum possible score value, independently of 
       % whether that is the worst of the best score.
       Smin = S.BasePoints;
       Smax = S.BasePoints;
       
       for i = 1 : NumPreds
          Smin = Smin + min(S.(PredNames{i}).Points);
          Smax = Smax + max(S.(PredNames{i}).Points);
       end
       
    end

    
end


end