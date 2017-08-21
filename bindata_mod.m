function bd = bindata_mod(obj,varargin)
% BINDATA Return a table of binned predictors.
%
% Syntax:
%
%      bdata = bindata(sc);
%      bdata = bindata(sc,data);
%      bdata = bindata(sc,data,param1,val1,...);
%
% Description:
%
%   BINDATA returns a table of binned predictors. Although BINDATA 
%   returns a table of the same size as the data input, only the predictors 
%   specified in the object's PredictorVars property are binned and the
%   remaining ones are unchanged.
%   
% Input Arguments:
%
%   sc - creditscorecard object.
%
% Optional Input Parameter Name/Value Pairs:
%
%   data - Data to bin, given the rules (cut points and category groupings)
%     set in the creditscorecard object. By default, it is set to the
%     object's raw data.
%
%   'OutputType'     Formats the output in the specified format:
%
%                    o 'BinNumber'     (Default) Returns the bin numbers 
%                                      corresponding to each observation.
%
%                    o 'Categorical'   Returns the bin label corresponding
%                                      to each observation.
%
%                    o 'WOE'           Returns the Weight Of Evidence
%                                      corresponding to each observation.
%
%                    o 'WOEModelInput' Use this option when fitting a
%                                      model.
%                                      o It returns the Weight Of Evidence
%                                        corresponding to each observation.
%                                      o It returns predictor variables,
%                                        but no IDVar or unused variables
%                                        are included in the output.
%                                      o It discards any predictors whose
%                                        bins have Inf or NaN WOE values.
%                                      o It includes the mapped response
%                                        variable as the last column.
%
%   'ResponseFormat' Specify whether the values of the raw response variable
%                    or the mapped response variable are used to bin the
%                    data.
% 
%                    o 'RawData' (Default) Uses the raw values of the
%                                response variable to bin the data.
%
%                    o 'Mapped'  Uses the mapped values of the response
%                                variable to bin the data. With this
%                                option, 'good' is mapped to 1, and 'bad'
%                                is mapped to 0. See also the 'GoodLabel'
%                                parameter in the CREDITSCORECARD
%                                constructor.
%
% Output:
%
%   bdata - Table of binned predictor variables. This is a table of the
%     same size as the data input (see exception in notes below), but only
%     the predictors specified in the creditscorecard object's
%     PredictorVars property are binned and the remaining ones are
%     unchanged.
%
% Notes:
%
%   o When 'OutputType' is set to 'WOEModelInput', bdata only contains
%     the columns corresponding to binned predictors whose bins do not have
%     Inf or NaN WOE values, and the mapped response as the last column.
%     Additionally, input data (when passed in) must contain the response 
%     variable when 'WOEModelInput' is used.
%
%   o Missing data (if any) are included in the output as missing data
%     as well, and do not influence the rules to discard predictors when
%     'WOEModelInput' is selected.
%
% See also BININFO, MODIFYBINS, AUTOBINNING, PLOTBINS.
%

% Copyright 2014-2015 The MathWorks, Inc.

    parser = inputParser();
            
    parser.addOptional('data',obj.RawData,@istable);

    parser.addParameter('OutputType','BinNumber',@(x)ismember(...
                lower(x),{'binnumber','categorical','woe','woemodelinput'}));

    parser.addParameter('ResponseFormat','RawData',@(x)ismember(...
                lower(x),{'rawdata','mapped'}));
            
    parser.parse(varargin{:});


    OutputType  = parser.Results.OutputType;
    ResponseFmt = parser.Results.ResponseFormat;
    bd          = obj.PreprocessPredictorData(parser.Results.data);

    hasResponse = ismember(obj.ResponseVar,bd.Properties.VariableNames);
    
    % If WOEModelInput, switch to mapped response
    if strcmpi(OutputType,'woemodelinput')
        
       if ~hasResponse
           error(message('finance:creditscorecard:bindata:MissingResponse'));
       end
        
       % Warn if ResponseFormat was provided AND it was set to RawData
       if ~ismember('ResponseFormat',parser.UsingDefaults)... 
             && strcmpi(ResponseFmt,'rawdata')
         
           warning(message('finance:creditscorecard:bindata:SwitchToMapped'));
       end
       
       ResponseFmt = 'mapped';

    end
    
    PredictorsInData = intersect(bd.Properties.VariableNames,...
       obj.PredictorVars,'stable');
   
    if isempty(PredictorsInData)
        error(message('finance:creditscorecard:bindata:InvalidPredictorVars'));
    end

    if hasResponse
        % Make sure the response variable in the input data takes
        % on the same values as the original data in the
        % creditscorecard object.
        ResponseVal = bd.(obj.ResponseVar);
        ResponseVal = unique(ResponseVal(~isnan(double(ResponseVal))));

        if ~all(ismember(ResponseVal,obj.ResponseOrder))
            error(message('finance:creditscorecard:bindata:InvalidResponseValues'));
        end
    end
    
    ModelInputVars = [PredictorsInData obj.ResponseVar];
    
    if strcmpi(OutputType,'woemodelinput')
       bd = bd(:,ModelInputVars);
    end

    PredictorsWithInf = cell(numel(obj.PredictorVars),1);
    
    for i = 1 : numel(PredictorsInData)
        
        PredName = PredictorsInData{i};
        
        [obj,Ind] = obj.getContainerIndex(PredName);
        
        if strcmpi(OutputType,'categorical')
            
            BinLabels = obj.BinContainers{Ind}.getLabels();
            bd.(PredName) = obj.BinContainers{Ind}.binData(...
                bd.(PredName),'Values',BinLabels);

        elseif strcmpi(OutputType,'woe') || strcmpi(OutputType,'woemodelinput')
            
            S = obj.BinContainers{Ind}.getStatistics('StatsList','woe');
            bd.(PredName) = obj.BinContainers{Ind}.binData(...
                bd.(PredName),'Values',S.woe);

            if any(~isfinite(S.woe))
                PredictorsWithInf{i} = PredName;
                if strcmpi(OutputType,'woemodelinput')
                   bd.(PredName) = [];
                end
            end
            
        else
            
            bd.(PredName) = obj.BinContainers{Ind}.binData(bd.(PredName));
            
        end

    end
       
    PredictorsWithInf(cellfun(@isempty,PredictorsWithInf)) = [];
    if ~isempty(PredictorsWithInf)
        
        if strcmpi(OutputType,'woemodelinput')
            warning(message('finance:creditscorecard:bindata:PredictorsWithInf1'));
            
        elseif strcmpi(OutputType,'woe')
            warning(message('finance:creditscorecard:bindata:PredictorsWithInf2'));
            
        end
        
        fprintf('\n')
        disp(PredictorsWithInf)
       
    end

    % Use the raw or mapped response variable values.
    if hasResponse
        if strcmpi(ResponseFmt,'mapped')
            ResponseValues = unique(bd.(obj.ResponseVar));
            [~,mappedRespVar] = obj.MapResponseVar(bd.(obj.ResponseVar),...
                                    ResponseValues);
            bd.(obj.ResponseVar) = mappedRespVar;

        end
    end

end
    
