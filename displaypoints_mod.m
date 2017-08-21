function [PointsInfo,MinScore,MaxScore] = displaypoints(obj)
% DISPLAYPOINTS Get creditscorecard points information per predictor, per bin.
%
% Syntax:
%
%   PointsInfo = displaypoints(sc);
%   [PointsInfo,MinScore,MaxScore] = displaypoints(sc);
%
% Description:
%
%   DISPLAYPOINTS returns a table of points for all bins of all predictor
%   variables used in the scorecard model. The table display information on
%   the predictor name, bin labels and the corresponding points per bin.
%   It also returns the minimum and maximum possible total scores.
%
% Input Arguments:
%
%   sc - creditcorecard object.
%
%
% Outputs:
%
%   PointsInfo - A table containing the creditscorecard points. That is, one
%      row for each predictor and each bin, with the corresponding points.
%
%         Predictors       Bin          Points
%         ----------       ---          ------
%
%         Predictor_1     Bin_11       Points_11
%                           ...           ...
%         Predictor_2     Bin_21       Points_21
%                           ...           ...
%         Predictor_j     Bin_ji       Points_ji
%                           ...           ...
%      When base points are reported separately (see FORMATPOINTS), the
%      first row of the table contains the base points.
%
%   MinScore - A scalar containing the minimum possible total score.
%
%   MaxScore - A scalar containing the maximum possible total score.
%
% Notes:
%
% o Minimum score is the lowest possible total score in the mathematical
%   sense, independently of whether a low score means high risk or low
%   risk. Similarly for maximum score.
%
% Algorithms:
%
%   The points for predictor j and bin i are, by default, given by
%
%      Points_ji = (Shift + Slope*b0)/p + Slope*(bj*WOEj(i))
%
%   where bj is the model coefficient of predictor j, p is the number of
%   predictors in the model, and WOEj(i) is the weight of evidence value
%   for the i-th bin corresponding to the j-th model predictor. 'Shift'
%   and 'Slope' are scaling constants.
%
%   The minimum and maximum scores are:
%
%      MinScore = Shift + Slope*b0 + min(Slope*b1*WOE1) + ... +min(Slope*bp*WOEp)),
%      MaxScore = Shift + Slope*b0 + max(Slope*b1*WOE1) + ... +max(Slope*bp*WOEp)).
%
%   Use FORMATPOINTS to control the way points are scaled, rounded, and
%   whether the base points are reported separately. See FORMATPOINTS
%   for more information on format parameters, and see the 'Algorithms'
%   section in SCORE for details and formulas on formatting options.
%
%   See also FORMATPOINTS, SCORE.

%   Copyright 2014 - 2015 The MathWorks, Inc.

    try
        [S,MinScore,MaxScore] = obj.getFormattedPoints;
    catch MException
        throw(MException)
    end
    
    BasePoints = S.BasePoints;
    S = rmfield(S,'BasePoints');
    S = rmfield(S,'ZeroWOEValue');

    Lengths = structfun(@(x)size(x,1),S);
    fields = fieldnames(S);

    Col1 = cell(sum(Lengths),1);
    Col2 = cell(sum(Lengths),1);
    Col3 = zeros(sum(Lengths),1);

    for i = 1 : numel(fields)
        if i == 1
            iStart = 1;
        else
            iStart = 1 + sum(Lengths(1:i-1));
        end
        iEnd   = sum(Lengths(1:i));    

        Predictors = repmat(S.(fields{i}).Properties.VariableNames(1),...
                     size(S.(fields{i}),1),1);
        Labels = S.(fields{i}).(fields{i});
        PointsInfo = S.(fields{i}).Points;
        Col1(iStart:iEnd) = Predictors;
        Col2(iStart:iEnd) = Labels;
        Col3(iStart:iEnd) = PointsInfo;

    end

    if obj.PointsFormat.BasePoints
       Col1 = ['BasePoints';Col1];
       Col2 = ['BasePoints';Col2];
       Col3 = [BasePoints;Col3];
    end
    
    PointsInfo = table(Col1,Col2,Col3,'VariableNames',{'Predictors',...
             'Bin','Points'});

    if strcmpi(obj.PointsFormat.Round,'finalscore')
       MinScore = round(MinScore);
       MaxScore = round(MaxScore);
    end

 end
   