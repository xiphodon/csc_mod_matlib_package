function obj = formatpoints_mod(obj,varargin)
% FORMATPOINTS Format creditscorecard points.
%
% Syntax:
%
%   sc = formatpoints(sc,param1,val1,...)
%
% Description:
%
%   Use FORMATPOINTS to modify the creditscorecard points: for example,
%   to change the scaling of the scores, or the rounding of the points. 
%   See the optional input parameters for a full list of point attributes 
%   that can be modified.
%
% Input Arguments:
%
%   sc - creditcorecard object.
%
% Optional Input Parameter Name/Value Pairs:
%
%   'BasePoints'         A logical scalar. If true, the scorecard explicitly
%                        separates base points. If false (default), the base 
%                        points are spread across all variables in the scorecard.
%                        See SCORE's 'Algorithms' section for details.
%
%   'Missing'            String to specify the points that should be
%                        assigned to missing or out-of-range information
%                        while scoring. Valid values are:
%
%                        o 'NoScore'    (Default) By default, missing and
%                                       out-of-range data do not get points
%                                       assigned, points are set to NaN.
%                                       The total score is also NaN.
%
%                        o 'ZeroWOE'    Missing or out-of-range data get
%                                       assigned a zero Weight-Of-Evidence
%                                       (WOE) value.
%
%                        o 'MinPoints'  Missing or out-of-range data get
%                                       the minimum possible points for
%                                       that predictor. This penalizes the
%                                       score if higher scores are better.
%
%                        o 'MaxPoints' Missing or out-of-range data get
%                                      the maximum possible points for
%                                      that predictor. This penalizes the
%                                      score if lower scores are better.
%
%   'Round'              String to indicate whether to round points or
%                        scores. Valid values are:
%
%                        o 'None'       (Default) No rounding is applied.
%
%                        o 'AllPoints'  Apply rounding to each predictor's 
%                                       points before adding up the total 
%                                       score.
%
%                        o 'FinalScore' Round the final score only (rounding 
%                                       is applied after all points have 
%                                       been added up).
%
% The Following Parameter Name/Value Pairs are Scaling Alternatives:
%
%   'PointsOddsAndPDO'   Numeric array with three elements [Points, Odds, PDO]
%                        that indicate target points (Points) for a given 
%                        odds level (Odds), and the desired number of points 
%                        to double the odds (PDO). Odds must be a positive 
%                        number. These values are used to find scaling 
%                        parameters for the scoring model. See SCORE's 
%                        'Algorithms' section for details.
%
%   'ShiftAndSlope'      Numeric array with two elements [Shift, Slope]
%                        that indicate the desired shift and slope scaling 
%                        parameters for the scorecard. Slope cannot be zero. 
%                        These values are used to scale the scoring model. 
%                        See SCORE's 'Algorithms' section for details.
%
%   'WorstAndBestScores' Numeric array with two elements [WorstScore, BestScore]
%                        that indicate the worst (highest risk) and best 
%                        (lowest risk) scores desired in the scorecard. 
%                        WorstScore and BestScore must be different. These 
%                        values are used to find scaling parameters for the
%                        scoring model. See SCORE's 'Algorithms' section for
%                        details.
%
% Output:
%
%   sc - An updated creditscorecard object.
%
%
% Notes:
%
%   o The scaling parameters 'PointsOddsAndPDO', 'ShiftAndSlope', and
%     'WorstAndBestScores' are mutually exclusive alternatives, hence only
%     one can be specified as optional input in a call to FORMATPOINTS.
%
%   o To remove a previous scaling and revert to unscaled scores, set
%     'ShiftAndSlope' to [0,1].
%
%   o 'Worst' score means the riskiest score, and its value could be lower or
%     higher than the 'best' score. In other words, the 'minimum' score may 
%     be the 'worst' score or the 'best' score, depending on the desired 
%     scoring scale.
%
%   o Similarly, the points to double the odds (PDO) may be positive or
%     negative, depending on whether higher scores mean lower risk, or vice
%     versa.
%
% See also SCORE, DISPLAYPOINTS.

% Copyright 2014 - 2015 The MathWorks, Inc.

    parser = inputParser();
    parser.addParameter('BasePoints',obj.PointsFormat.BasePoints,@islogical);

    parser.addParameter('Round',obj.PointsFormat.Round,@(x)ismember(lower(x),...
                        {'none','allpoints','finalscore'}));

    parser.addParameter('Missing',obj.PointsFormat.Missing,@(x)ismember(lower(x),...
       {'noscore','zerowoe','zeropoints','minpoints','maxpoints'}));

    parser.addParameter('ShiftAndSlope',[],...
                @(x)(isempty(x)||(isnumeric(x) && numel(x)==2)));

    parser.addParameter('WorstAndBestScores',[],...
                @(x)(isempty(x)||(isnumeric(x) && numel(x)==2)));

    parser.addParameter('PointsOddsAndPDO',[],...
                @(x)(isempty(x)||(isnumeric(x) && numel(x)==3)));

    parser.parse(varargin{:});

    obj.PointsFormat.BasePoints = parser.Results.BasePoints;
    obj.PointsFormat.Round      = parser.Results.Round;

    Missing = parser.Results.Missing;
    % Reset 'ZeroPoints' option in R2014b to 'ZeroWOE'
    if strcmpi(Missing,'zeropoints')
        Missing = 'ZeroWOE';
    end
    obj.PointsFormat.Missing = Missing;
    
    ShiftAndSlope = parser.Results.ShiftAndSlope;
    WorstAndBestScores = parser.Results.WorstAndBestScores;
    PointsOddsAndPDO = parser.Results.PointsOddsAndPDO;
    
    if ~isempty(WorstAndBestScores)
       
       % Only one scaling method can be used
       if ~isempty(ShiftAndSlope) || ~isempty(PointsOddsAndPDO)
          error(message('finance:creditscorecard:formatpoints:MultipleScalingInputs'))
       end
       
       if WorstAndBestScores(1) == WorstAndBestScores(2)
            error(message('finance:creditscorecard:formatpoints:WorstBestAreEqual'))
       end
       
       obj.PointsFormat.ShiftAndSlope = [];
       obj.PointsFormat.WorstAndBestScores = WorstAndBestScores;
       obj.PointsFormat.PointsOddsAndPDO = [];

    elseif ~isempty(PointsOddsAndPDO)
       
       % Only one scaling method can be used
       % If we enter this 'if', we know WorstAndBestScores is empty, but
       % still need to check that ShiftAndSlope was not provided
       if ~isempty(ShiftAndSlope)
          error(message('finance:creditscorecard:formatpoints:MultipleScalingInputs'))
       end
       
       if PointsOddsAndPDO(2)<=0 % Odds must be positive
          error(message('finance:creditscorecard:formatpoints:InvalidOdds'));
       end
       
       obj.PointsFormat.ShiftAndSlope = [];
       obj.PointsFormat.WorstAndBestScores = [];
       obj.PointsFormat.PointsOddsAndPDO = PointsOddsAndPDO;

    elseif ~isempty(ShiftAndSlope)
       
       % If we enter here, we know PointsOddsAndPDO and WorstAndBestScores
       % were empty, no need to validate that only ShiftAndSlope was given
       
       if ShiftAndSlope(2)==0 % Slope cannot be zero
          error(message('finance:creditscorecard:formatpoints:InvalidSlope'));
       end
       
       obj.PointsFormat.ShiftAndSlope = ShiftAndSlope;
       obj.PointsFormat.WorstAndBestScores = [];
       obj.PointsFormat.PointsOddsAndPDO = [];

    end

end