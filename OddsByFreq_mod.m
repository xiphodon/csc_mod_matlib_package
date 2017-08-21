function odds = OddsByFreq_mod(FreqTab, NanFreqTab)
% ODDSBYFREQUENCY Compute odds given a frequency table.
%
%   This is an internal function that is not meant to be called directly
%   by the user.

% Copyright 2014 The MathWorks, Inc. 

if size(FreqTab,2)~=2
   error('Only binary tables (two columns) are supported');
end

odds = FreqTab(:,1)./FreqTab(:,2);
if ~isempty(NanFreqTab) && any(NanFreqTab ~= 0)
    odds(end + 1) = NanFreqTab(1) / NanFreqTab(2);
end