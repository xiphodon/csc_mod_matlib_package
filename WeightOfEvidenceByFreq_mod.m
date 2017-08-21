function WOE = WeightOfEvidenceByFreq_mod(FreqTab, NanFreqTab)
% WEIGHTOFEVIDENCEBYFREQUENCY Compute Weight Of Evidence (WOE) given a
%   frequency table.
%
%   This is an internal function that is not meant to be called directly
%   by the user.

% Copyright 2014 The MathWorks, Inc. 

if size(FreqTab,2)~=2
   error('Only binary tables (two columns) are supported');
end

nPerCol = sum(FreqTab, 1);
if ~isempty(NanFreqTab)
    nPerCol = nPerCol + NanFreqTab;
end
p = bsxfun(@rdivide,FreqTab,nPerCol);

WOE = log(p(:,1)./p(:,2));
if ~isempty(NanFreqTab) && any(NanFreqTab ~= 0)
    WOE(end + 1) = log((NanFreqTab(1)/ nPerCol(1))/(NanFreqTab(2)/nPerCol(2)));
end
