classdef EqualFrequency_mod < Algorithm_mod
% EQUALFREQUENCY Equal-frequency automatic binning algorithm.
%
% Description:
%
%   This is an unsupervised automatic binning algorithm that divides the
%   data into a predetermined number of bins that contain approximately the
%   same number of observations. This algorithm is also known as "equal
%   height" or "equal depth".
%
%   Let v[1], v[2],..., v[N] be the sorted list of different values or
%   categories observed in the data. Let f[i] be the frequency of v[i]. Let
%   F[k] = f[1]+...+f[k] be the cumulative sum of frequencies up to the
%   k-th sorted value. Then F[N] is the same as the total number of
%   observations.
%
%   Define AvgFreq = F[N] / NumBins, which is the ideal or target average
%   frequency per bin. Similarly, the target cumulative frequency for the
%   first n bins is n*AvgFreq. The implementation of the equal-frequency
%   algorithm attempts to match the cumulative frequency up to the n-th
%   bin, and sets the n-th cut point by finding the index k such that the
%   distance abs(F[k] - n*AvgFreq) is minimized.
%
%   If a single value contains too many observations, equal-frequency bins
%   are not possible, and the above rule yields less than NumBins total
%   bins. In that case, the algorithm determines NumBins bins by breaking
%   up bins, in the order in which the bins were constructed.
%
%   This is an internal function that is not meant to be called directly
%   by the user.

% Copyright 2014 The MathWorks, Inc.

   properties (Access = protected)
   end

   methods
      
      function obj = EqualFrequency_mod(algoID)
         % EQUALFREQUENCY Create an equal-frequency algorithm object.
         %
         % The algorithm ID is an optional input.
         
         if nargin<1
            algoID = 'AlgoEqFreq';
         end
         obj = obj@Algorithm_mod(algoID);
         
         % This is an unsupervised algorithm
         obj.IsSupervised = false;
         
         % Supported algorithm options
         obj = obj.addOption('NumBins',5,@(x)(x>=1&mod(x,1)==0));
	
         obj = obj.addOption('SortCategories','odds',@(x)ismember(lower(x),...
                             {'odds','goods','bads','totals','none'}));
         
      end
      
      function dbc = runAlgorithm(obj,dbc)
         % RUNALGORITHM Apply equal-frequency binning to given data.
         %
         % Inputs:
         %
         %   o obj: EqualFrequency algorithm object.
         %
         %   o dbc: Data binning container object. This must be either a
         %          Numeric or Categorical object derived from
         %
         %           internal.finance.binning.container.Container
         %

         if ~isa(dbc,'Container_mod')
            error(message('finance:internal:finance:binning:EqualFrequency:InvalidContainer'));
         end
         
         NumBins = obj.getOptions('NumBins');
         
         % Pre-processing of the data: apply sorting for categorical 
         % nominal predictors only
         if isa(dbc,'Categorical_mod')
             if ~dbc.isOrdinal
                 Sorting = obj.getOptions('SortCategories');
                 dbc = obj.applySorting(dbc,Sorting);
             end
         end
         
         FreqTab = dbc.getFrequencyTable('BinLevel',false);
         Freq    = sum(FreqTab,2);
         CumFreq = cumsum(Freq);
         NumObs  = CumFreq(end);
         NumRows = size(FreqTab,1);

         %add by zhenxiao
         if NumRows <= NumBins
            % Not enough data values / levels to split into desired number
            % of bins; keep each data value / level in its own bin
            CutPointIndices = 2:NumRows;
         else
            % This approach sets a "cumulative target" and fills in
            % independently of what happened in previous bins. The cut
            % point for the i+1 bin is set at the value where the
            % cumulative frequency is closer to the cumulative target.
            
                AvgObsPerBin = NumObs / NumBins;
                CutPointIndices = zeros(1,NumBins-1);
                for ii=1:NumBins-1
                   CumTarget = ii*AvgObsPerBin;
                   ind = find(CumFreq>CumTarget,1,'first');
                   if ind == 1
                      CutPointIndices(ii) = 2;
                   else
                      if abs(CumFreq(ind-1)-CumTarget) < abs(CumFreq(ind)-CumTarget)
                         CutPointIndices(ii) = ind;
                      else
                         if ind < NumRows
                            CutPointIndices(ii) = ind+1;
                         else
                            CutPointIndices(ii) = ind;
                         end
                      end
                   end
                end
            % Handling of small number of rows. It can lead to repeated
            % cutpointindices, or cutpointindex equal to 1. If duplicates
            % are eliminated, fill in with individual rows in the order
            % they appear.
            CutPointIndices = unique(max(2,CutPointIndices));
            % modified by zhenxiao
            %if length(CutPointIndices)<NumBins-1
            %   n = (NumBins-1) - length(CutPointIndices);
            %   UnusedInd = setdiff(2:NumRows,CutPointIndices);
            %   CutPointIndices = sort([CutPointIndices UnusedInd(1:n)]);
            %end
         end
         
         dbc = obj.setBinsByCutPointIndices(dbc,CutPointIndices, 0);

      end
      
   end
   
end
