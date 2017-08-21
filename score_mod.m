function [Scores,Points] = score_mod(obj,varargin)
% SCORE Determine credit score for given data.
%
% Syntax:
%
%   Scores = score(sc)
%   [Scores,Points] = score(sc,data)
%
% Description:
%
%   SCORE computes the credit scores for the given data. If the data
%   input is not explicitly provided, SCORE determines scores for the
%   CREDITSCORECARD data.
%
% Optional Input Parameter Name/Value Pairs:
%
%   data - Table of new data (e.g test data). Each row corresponds to
%       individual observations. The data must contain columns for each of 
%       the predictors in the credit scorecard model.
%
% Output:
%
%  Scores  -  Vector of scores for each observation.
%
%  Points  -  Table of points per predictor for each observation.
%
% Notes:
%
% o FORMATPOINTS supports multiple alternatives to modify the scaling
%   of the scores.
%
% o FORMATPOINTS can also be used to control the rounding of points and
%   scores, and whether the base points are to be reported separately.
%
% o Missing data translates into NaN values for the corresponding
%   points, and therefore for the total score. Use FORMATPOINTS to
%   modify the score behavior for rows with missing data.
%
%
% Algorithms:
%
%   The score of an individual i is given by the formula
%
%      Score(i) = Shift + Slope*(b0 + b1*WOE1(i) + b2*WOE2(i)+ ... +bp*WOEp(i))
%
%   where bj is the coefficient of the j-th variable in the model, and
%   WOEj(i) is the Weight Of Evidence (WOE) value for the i-th individual
%   corresponding to the j-th model variable. Shift and Slope are scaling
%   constants further discussed below. The scaling constant can be
%   controlled with FORMATPOINTS.
%
%   If the data for individual i is in the i-th row of a given dataset, to
%   compute a score, the data(i,j) is binned using existing binning maps,
%   and converted into a corresponding Weight Of Evidence value WOEj(i).
%   Using the model coefficients the _unscaled_ score is computed as
%
%      s = b0 + b1*WOE1(i) + ... +bp*WOEp(i).
%
%   For simplicity, we assume in the description above that the j-th
%   variable in the model is the j-th column in the data input, although,
%   in general, the order of variables in a given dataset does not have to
%   match the order of variables in the model, and the dataset could have
%   additional variables that are not used in the model.
%
%   The formatting options can be controlled using the FORMATPOINTS
%   method. See FORMATPOINTS for a description of the syntax to modify the 
%   formatting of points. The remainder of this section contains a mathematical 
%   description of how the formatting options work.
%
%   When the base points are reported separately (see FORMATPOINTS
%   parameter 'BasePoints'), the base points are given by
%
%      Base Points = Shift + Slope*b0,
%
%   and the points for the j-th predictor, i-th row are given by
%
%      Points_ji = Slope*(bj*WOEj(i))).
%
%   By default, the base points are not reported separately, in which
%   case
%
%      Points_ji = (Shift + Slope*b0)/p + Slope*(bj*WOEj(i)),
%
%   where p is the number of predictors in the creditscorecard model.
%
%   By default, no rounding is applied to the points ('Round' is 'None').
%   If 'Round' is set to 'AllPoints', then the points for individual i for
%   variable j are given by
%
%      points if rounding is 'AllPoints':round( Points_ji )
%
%   and, if base points are reported separately, the are also rounded.
%   This yields integer-valued points per predictor, hence also integer-
%   valued scores. If 'Round' is set to 'FinalScore', then the points per 
%   predictor are not rounded, and only the final score is rounded
%
%      score if rounding is 'FinalScore': round( Score(i) ).
%
%   Regarding the scaling parameters, the Shift parameter, and the Slope 
%   parameter can be set directly with the 'ShiftAndSlope' parameter of 
%   FORMATPOINTS.
%
%   Alternatively, the user can use FORMATPOINTS' 'WorstAndBestScores'
%   parameter. In this case, the parameters Shift and Slope are found
%   internally by solving the system
%
%      Shift + Slope*smin = WorstScore,
%      Shift + Slope*smax = BestScore,
%
%   where WorstScore and BestScore are the first and second elements in
%   FORMATPOINTS' 'WorstAndBestScores' parameter, and smin and smax are
%   the minimum and maximum possible unscaled scores:
%
%      smin = b0 + min(b1*WOE1) + ... +min(bp*WOEp),
%      smax = b0 + max(b1*WOE1) + ... +max(bp*WOEp).
%
%   A third alternative to scale scores is the 'PointsOddsAndPDO'
%   parameter in FORMATPOINTS. In this case, it is assumed that the
%   unscaled score s gives the log-odds for a row, and the Shift and
%   Slope parameters are found by solving the following system
%
%      Points = Shift + Slope*log(Odds)
%      Points + PDO = Shift + Slope*log(2*Odds)
%
%   where Points, Odds, and PDO ("points to double the odds") are the
%   first, second, and third elements in the 'PointsOddsAndPDO'
%   parameter.
%
%   Last, whenever a given dataset has a missing or out-of-range value
%   data(i,j), the points for predictor j, for individual i, are set to NaN
%   by default, which results in a missing score for that row (a NaN
%   score). Using FORMATPOINT's 'Missing' parameter, the user can modify
%   this behavior and set the corresponding Weight-Of-Evidence (WOE) value
%   to zero, or set the points to the minimum points or the maximum points
%   for that predictor.
%
%   See also FORMATPOINTS, DISPLAYPOINTS.

% Copyright 2014-2015 The MathWorks, Inc.

    parser = inputParser();
    parser.addOptional('data',obj.RawData,@istable);

    parser.parse(varargin{:});
    
    data = obj.PreprocessPredictorData(parser.Results.data);

    if isempty(obj.BinContainers)
        error(message('finance:creditscorecard:score:NonBinnedData'));
    end
    
    if isempty(obj.ModelVars) || isempty(obj.ModelCoeff)
        error(message('finance:creditscorecard:score:EmptyModel'));
    end
    
    if (numel(obj.ModelCoeff) == 1 && strcmpi(obj.ModelVars{:},'(intercept)'))
        error(message('finance:creditscorecard:score:InterceptOnlyModel'));
    end
    
    S = obj.getFormattedPoints;
    
    FieldNames = fields(S);
    PredNames = setdiff(FieldNames,{'BasePoints','ZeroWOEValue'},'stable');
    NumPreds = length(PredNames);
    
    if any(~ismember(PredNames,data.Properties.VariableNames))
       error(message('finance:creditscorecard:score:MissingPredictor'));
    end
    
    Points = data(:,PredNames);
    
    for i = 1:NumPreds
       
       PointsPred = S.(PredNames{i}).Points;
       
       switch lower(obj.PointsFormat.Missing)
          % case 'noscore': no need to change anything in this case
          case 'zerowoe'
             PointsPred = [PointsPred(:); S.ZeroWOEValue];
             
          case 'minpoints'
             PointsPred = [PointsPred(:); min(PointsPred)];
             
          case 'maxpoints'
             PointsPred = [PointsPred(:); max(PointsPred)];
             
       end
       
       [obj,Ind] = obj.getContainerIndex(PredNames{i});
       
       Points.(PredNames{i}) = obj.BinContainers{Ind}.binData(data.(PredNames{i}),...
          'Values',PointsPred);
       
    end
    
    if obj.PointsFormat.BasePoints
       Points.BasePoints = S.BasePoints*ones(size(data,1),1);
       NumCols = size(Points,2);
       Points = Points(:,[NumCols 1:NumCols-1]);
    end
    
    Scores = zeros(size(data,1),1);
    
    if obj.PointsFormat.BasePoints
       Scores = Scores + Points.BasePoints;
    end
    
    for i = 1:NumPreds
       Scores = Scores + Points.(PredNames{i});
    end
    
    if strcmpi(obj.PointsFormat.Round,'finalscore')
       Scores = round(Scores);
    end

end
    
