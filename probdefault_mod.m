function pd = probdefault_mod(obj,varargin)
% PROBDEFAULT computes the probability of default.
%
% Syntax:
%
%      pd = probdefault(sc);
%      pd = probdefault(sc,data);
%
% Description:
%
%   PROBDEFAULT computes the likelihood of default for a given dataset. By
%   default, the data used to build the credit scorecard object are used. 
%   The user can also supply input data, to which the same computation of 
%   probablility of default will be applied.
%
%   
% Input Arguments:
%
%   sc - creditscorecard object.
%
% Optional Input Parameter Name/Value Pairs:
%
%   data - Data to which the same probability of default rules will be 
%     applied. By default, this argument is set to the object's training
%     data.
%
% Output:
%
%   pd - NumObs-by-1 numerical array of default probabilities.
%
%
% Notes:
%
%   Once the unscaled points are computed, the probability of being "Good" 
%   is given by the following formula (logistic function):
%
%       ProbGood = 1./(1 + exp(-UnscaledScores)).
% 
%   The probability of default is thus:
%
%       pd = 1 - ProbGood;
%
%   For the definition of unscaled scores, see the functions FORMATPOINTS, 
%   SCORE.
%
%   See also FORMATPOINTS, SCORE.


% Copyright 2015 The MathWorks, Inc.

    parser = inputParser();
    
    parser.addOptional('data',obj.RawData,@istable);
    
    parser.parse(varargin{:});
    
    data = parser.Results.data;

    % Format points to get the unscaled points.
    obj = formatpoints_mod(obj, 'ShiftAndSlope',[0 1],'Round','None');

    UnscaledScores = score_mod(obj,data);

    % Because 'Good' is mapped to '1' during the logistic fitting, the logit
    % transformation yields probabilities of being good. The probability of
    % default is thus 1 minus that probability.

    ProbGood = 1./(1 + exp(-UnscaledScores));

    pd  = 1 - ProbGood;

end
