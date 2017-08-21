function Stats = GetStatisticsByFreq_mod(FreqTab,NanFreqTab,varargin)
% GETSTATISTICSBYFREQUENCY Compute statistics given a frequency table.
%
%   This is an internal function that is not meant to be called directly
%   by the user.

% Copyright 2014 - 2015 The MathWorks, Inc. 

parser = inputParser;
parser.addParameter('StatsList',{'all'},@(x)all(ismember(lower(x),...
   {'all','entropy','infoval','odds','woe'})));

parser.parse(varargin{:});

StatsList = cellfun(@(x)lower(x),parser.Results.StatsList,...
   'UniformOutput',false);

if any(ismember({'all','entropy'},StatsList))
   [e,erows] = internal.finance.binning.utils.EntropyPartitionByFreq(FreqTab);
   Stats.entropy = e;
   Stats.entropybins = erows;
end

if any(ismember({'all','infoval'},StatsList))
   [iv,ivrows] = InformationValueByFreq_mod(FreqTab, NanFreqTab);
   Stats.infoval = iv;
   Stats.infovalbins = ivrows;
end

if any(ismember({'all','odds'},StatsList))
   Stats.odds = OddsByFreq_mod(FreqTab, NanFreqTab);
end

if any(ismember({'all','woe'},StatsList))
   Stats.woe = WeightOfEvidenceByFreq_mod(FreqTab, NanFreqTab);
end

