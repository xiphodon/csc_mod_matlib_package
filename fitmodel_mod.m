function [obj,mdl] = fitmodel(obj,varargin)
% FITMODEL Fit logistic regression model to Weight Of Evidence (WOE) data.
%
% Syntax:
%
%    [sc,mdl] = fitmodel(sc) 
%    [sc,mdl] = fitmodel(sc,param1,val1,...)
%
% Description:
%
%   FITMODEL fits a linear logistic regression model to the Weight Of
%   Evidence (WOE) data, and stores the model predictor names and 
%   corresponding coefficients in the creditscorecard object.
%
% Input Arguments:
%
%   sc - creditcorecard object.
% 
% Optional Input Parameter Name/Value Pairs:
%
%   'PredictorVars'     Cell array of strings to be used as predictors in
%                       FITMODEL. When this parameter is provided, the
%                       'PredictorsVars' property of the creditscorecard
%                       object is updated. If not provided, the
%                       'PredictorVars' property of the creditscorecard
%                       object is used.
%
%   'VariableSelection' String that speficies the variable selection method
%                       to fit the logistic regression model. Options are:
%
%                       o 'Stepwise'  (Default) Use a stepwise procedure to
%                                     fit a linear regression model. Only
%                                     variables in 'PredictorVars' can
%                                     potentially become part of the model.
%                                     Calls 'stepwiseglm'. Use the
%                                     'StartingModel' parameter to select
%                                     the starting model.
%
%                       o 'FullModel': Fits a model with all predictor
%                                      variables in 'PredictorVars'. Calls
%                                      'fitglm'.
%
%   'StartingModel'      String that specifies the starting model for the
%                        'Stepwise' variable selection method. Options are: 
%
%                        o 'Constant' (Default) Start the stepwise method
%                                     with an empty (constant only) model.
%
%                        o 'Linear'   Start the stepwise method from a full
%                                     (all predictors in) model.
%
%                        This parameter is used only for the 'Stepwise'
%                        option of 'VariableSelection' and has no effect
%                        for the 'FullModel' option.
%
%   'Display'            String to indicates whether the model information
%                        is displayed at the command line. Options are:
%
%                        o 'On'  (Default) Display model information.
%
%                        o 'Off' Do not display model information.
%
% Ouput:
%
%   sc - The updated creditscorecard object, which contains information
%      about the model predictors and coefficients of the fitted model.
%
%   mdl - Object of type GeneralizedLinearModel containing the fitted
%      model.
%
% Notes:
%
%   o FITMODEL internally calls BINDATA to apply the Weight Of Evidence
%     transformation of the predictor variables. It then fits a linear
%     logistic regression model.
%
%   o The response variable is mapped so that 'good' is 1 and 'bad' is 0.
%
%   o Only linear terms are included in the model, no interactions or any
%     other higher-order terms.
%
%   o SETMODEL offers an alternative to FITMODEL. Users can get the Weight
%     Of Evidence data using BINDATA (use the 'WOEModelInput' option for
%     the 'OutputType' parameter in BINDATA), fit a linear logistic
%     regression model, and provide the names of the predictors in the
%     final model, along with their coefficients using SETMODEL. See
%     SETMODEL for more details.
%
% See also SETMODEL, STEPWISEGLM, FITGLM.

% Copyright 2014 The MathWorks, Inc.

    parser = inputParser();
    parser.addParameter('PredictorVars',obj.PredictorVars,@(x)iscell(x));
    
    parser.addParameter('VariableSelection','stepwise',...
                @(x)ismember(lower(x),{'stepwise','fullmodel'}));

    parser.addParameter('StartingModel','constant',...
                @(x)ismember(lower(x),{'constant','linear'}));   

    parser.addParameter('Display','On',@(x)ismember(lower(x),{'on','off'}));

    parser.parse(varargin{:});
    
    obj.PredictorVars = obj.PreprocessPredictorVar(parser.Results.PredictorVars);
                                
    VarSelection      = parser.Results.VariableSelection;
    StartingModel     = parser.Results.StartingModel;
    Display           = parser.Results.Display;

    woedata = bindata_mod(obj, 'OutputType','WOEModelInput');
    
    if strcmpi(VarSelection,'fullmodel')
        mdl = fitglm(woedata,'linear','Distribution','binomial');
        
    else
        %zhenxiao
        mdl = stepwiseglm(woedata,StartingModel,'Distribution',...
            'binomial','Upper','linear', 'PEnter',0.025,'PRemove',0.05);
        
    end
    
    obj.ModelVars  = mdl.CoefficientNames;
    obj.ModelCoeff = mdl.Coefficients.Estimate';
    
    if strcmpi(Display,'on')
        disp(mdl)
    end
    
end
    
   