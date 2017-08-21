function [InfoVal,InfoValRows] = InformationValueByFreq_mod(FreqTab, NanFreqTab)
% INFORMATIONVALUEBYFREQUENCY Compute information value given a frequency
%   table.
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
if ~isempty(NanFreqTab) && any(NanFreqTab ~= 0)
    nanp = NanFreqTab ./ nPerCol;
end

InfoValRows = (p(:,1)-p(:,2)).*log(p(:,1)./p(:,2));
InfoVal = nansum(InfoValRows);
if ~isempty(NanFreqTab) && any(NanFreqTab ~= 0)
    NanInfoValRow = (nanp(:,1)-nanp(:,2)).*log(nanp(:,1)./nanp(:,2));
    if ~isnan(NanInfoValRow)
        InfoValRows(end+1) = NanInfoValRow;
        InfoVal = InfoVal + NanInfoValRow;
    end
end

