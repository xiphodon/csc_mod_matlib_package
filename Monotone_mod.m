classdef Monotone_mod < Algorithm_mod
% MONOTONE Performs Monotone Adjacent Pooling algorithm (MAPA) to bin data.
%
% Description:
%
%   This is a supervised automatic binning algorithm, belonging to the
%   family of Pooling algorithms, used for coarse classing data. The 
%   Monotone Adjacent Pooling Algorithm (MAPA) assumes that there is come
%   sort of order in the data, and hence, that neighboring bins can be
%   grouped (merged) together. MAPA can thus be used on numeric data 
%   (discrete or continuous), as well as categorical data.
%
%   For either numeric and categorical data, the first step in binning 
%   involves:
%
%       - Computing the frequency table, FreqTab, of the data. FreqTab has
%         two columns corresponding to 'Good' and 'Bad' observations, 
%         respectively.
%       - Computing the total frequencies Freq defines as the sum, across
%         columns of FreqTab.
%
%   For numeric predictors, no pre-processing is needed and the MAPA works
%   as follows:
%       - Set the counter index iStart = 1.
%
%       - While iStart < NumBins, where NumBins is the lesser of the 
%         starting number of bins for the algorithm and of the number of 
%         observations for the given predictor:
%
%               * Compute the cumulative Good ratio, CumGoodRatio
%
%                     CumGoodRatio = cumsum(FreqTab(:,1))./cumsum(Freq)  
%
%                 and find the index, Incr, where either the minimum or 
%                 maximum (see value for 'Trend') occurs.
%
%               * Set the new value of iStart = iStart+Incr and define and
%                 and store the index of the cut-point iCutPoint =
%                 iStart-1.
%
%       - Once all cut-points are defined, set them to the binning
%         container object and return this latter.
%
%   For categorical predictors, since no intuitive order exists, some
%   pre-processing is needed:
%
%       - If some categories have zero entries for either 'Good' or 'Bad', 
%         they lead to Inf values for the WOE. In order to avoid initial WOE
%         curves with discontinuities, we first group together the
%         categories with zero count of 'Good' and those with zero count of
%         'Bad'. Then, we set their bin number (BinNum) as the last in the
%         CatGrouping property of the binning container.
%
%       - Sort the categories according to the criterion set in the
%         'Sorting' option (see above).
%   
%       - Follow the same steps as for numeric predictors above.
%
%   This is an internal function that is not meant to be called directly
%   by the user.
%

% Copyright 2015 The MathWorks, Inc.

   properties (Access = protected)
   end

   methods
      
      function obj = Monotone_mod(algoID)
         % MONOTONE Create a Monotone Adjacent Pooling Algorithm object.
         %
         % The algorithm ID is an optional input.
         
         if nargin<1
            algoID = 'AlgoMonotone';
         end
         obj = obj@Algorithm_mod(algoID);
         
         % This is a supervised algorithm
         obj.IsSupervised = true;
         
         % Supported algorithm options:
         % 'InitialNumBins': Initial number of bins that the algorithm
         %     starts with to bin data. Pre-processing applies Equal
         %     Frequency with InitialNumBins bins. The MAPA is applied. This
         %     algorithm option is relevant for numeric containers only. 
         %
         % 'Trend': This option sets the monotonicity trend of the Weight
         %     of Evidence (WOE) curve. The values are 'Increasing' and
         %     'Decreasing'.
         %
         % 'SortCategories': Sort categories, in ascending order of the
         %     selected criterion: 'Odds', 'Goods', 'Bads', 'Totals'. If no
         %     sorting is required, the value can be set to 'None'. This 
         %     algorithm option is used for categorical containers only.
         
         obj = obj.addOption('InitialNumBins',10,@(x)(x>=0&mod(x,1)==0));
         
         obj = obj.addOption('Trend','auto',@(x)(ischar(x)&&...
                             ismember(lower(x),{'auto','increasing','decreasing'})));
                         
         obj = obj.addOption('SortCategories','odds',@(x)ismember(lower(x),...
                             {'odds','goods','bads','totals','none'}));
         
      end
      
      
      function dbc = runAlgorithm(obj,dbc)
         % RUNALGORITHM Apply binning by the pooling method to given data.
         %
         % Inputs:
         %
         %   o obj: Monotone algorithm object.
         %
         %   o dbc: Data binning container object. This must be either a
         %          Numeric or Categorical object derived from
         %
         %           internal.finance.binning.container.Container
         %

         if ~isa(dbc,'Container_mod')
            error(message('finance:internal:finance:binning:Monotone:InvalidContainer'));
         end

         % If the initial number of bins is passed, it needs to be greater
         % than two.
         if obj.getOptions.InitialNumBins <= 2
             error(message('finance:internal:finance:binning:Monotone:InvalidInitNumBins'));
         end
         
         if isa(dbc,'Numeric_mod')
             dbc = obj.runMonotoneAdjacentPoolingNumeric(dbc);
             
         elseif isa(dbc,'Categorical_mod')
             dbc = obj.runMonotoneAdjacentPoolingCategorical(dbc);
         end
         
      end
      
   end
   
   
   methods ( Access = protected )
       
      function dbc = runMonotoneAdjacentPoolingNumeric(obj,dbc)
          % This method applies the Monotone Adjacent Pooling Algorithm 
          % (MAPA) to numeric variables. The pre-processing includes
          % applying an Equal Frequency algorithm first.
          
          % Step 1: Pre-processing. Create and apply an Equal Frequency
          % algorithm on the binning container object.
          % Note that in the event that the frequency table has less than
          % three rows, then EqualFrequency treats the cases where the
          % desired number of bins (InitialNumBins) is less than
          % modified by zhenxiao
          % InitNumBins = obj.getOptions('InitialNumBins');
    
          Trend = lower(obj.getOptions('Trend'));
          max_iv = 0.0;
          max_InitNumBins = 3;
          max_trend_idx = 1;
          Trends = {'increasing', 'decreasing'};
          for trend_idx = 1:2
              Trend = Trends{trend_idx};
          for InitNumBins = 3:20

          EqFreqAlgo = EqualFrequency_mod();
          EqFreqAlgo = EqFreqAlgo.setOptions('NumBins',InitNumBins);

          % Step 2: Pre-processing. Run the Equal Frequency algorithm on 
          % the binning container object
          dbc = EqFreqAlgo.runAlgorithm(dbc);

          % Step 3: Set up the frequency table that should be passed to the
          % main Monotone algorithm. For Numeric, we need the frequency
          % table at the bin level, since Equal Frequency was applied
          % above. 
          [FreqTab, NanFreqTab] = dbc.getFrequencyTable();
          % Step 4: Compute cut points then update the binning container 
          % 'dbc'.
              
          [CutPointIndices, nanBin] = obj.runMonotoneAdjacentPoolingMainForNumeric(FreqTab,NanFreqTab,...
              Trend);
          
          % Step 5: If 'CutPoints' is empty, do nothing. This will return
          % the container object with the Equal Frequency binning only.
          % if isempty(CutPointIndices) && nanBin ~= 0
              % continue
          % end
          
          %if ~isempty(CutPointIndices) || nanBin == 0
              EqFreqCutPoints = dbc.getCutPoints;
              CutPoints = EqFreqCutPoints(CutPointIndices-1);
              dbc = dbc.setCutPoints(CutPoints, nanBin);
          %end
          
          [new_FreqTab, new_NanFreqTab] = dbc.getFrequencyTable();
          [iv,ivrows] = InformationValueByFreq_mod(new_FreqTab, new_NanFreqTab);
          
          assert(~isinf(iv));  

          if iv > max_iv
            max_iv = iv;
            max_InitNumBins = InitNumBins;
            max_trend_idx = trend_idx;
          end
          %disp(['bin ' num2str(InitNumBins) ' iv ' num2str(iv)  ' mav_iv ' num2str(max_iv)])
          end
          end
          Trend = Trends{max_trend_idx};
          %disp(['max_bin ' num2str(max_InitNumBins) ' mav_iv ' num2str(max_iv)])
          EqFreqAlgo = EqualFrequency_mod();
          EqFreqAlgo = EqFreqAlgo.setOptions('NumBins',max_InitNumBins);

          % Step 2: Pre-processing. Run the Equal Frequency algorithm on 
          % the binning container object
          dbc = EqFreqAlgo.runAlgorithm(dbc);

          % Step 3: Set up the frequency table that should be passed to the
          % main Monotone algorithm. For Numeric, we need the frequency
          % table at the bin level, since Equal Frequency was applied
          % above. 
          [FreqTab, NanFreqTab] = dbc.getFrequencyTable();

          % Step 4: Compute cut points then update the binning container 
          % 'dbc'.
                           
          [CutPointIndices, nanBin] = obj.runMonotoneAdjacentPoolingMainForNumeric(FreqTab,NanFreqTab,...
              Trend);
          
          % Step 5: If 'CutPoints' is empty, do nothing. This will return
          % the container object with the Equal Frequency binning only.
          %if ~isempty(CutPointIndices) || nanBin == 0
              EqFreqCutPoints = dbc.getCutPoints;
              CutPoints = EqFreqCutPoints(CutPointIndices-1);
              dbc = dbc.setCutPoints(CutPoints, nanBin);
          %end
      end
      
      
      function dbc = runMonotoneAdjacentPoolingCategorical(obj,dbc)
          % This method applies the Monotone Adjacent Pooling Algorithm 
          % (MAPA). For categorical predictors, the bins are first "sorted" 
          % according to a set of criteria. Then, the MAPA is applied. In 
          % case the bins have Inf in the WOE data, then these bins are 
          % initially grouped together in the last bin.
          
          % Step 1: Set up the frequency table that should be passed to the
          % main Monotone algorithm. For categorical, we need the frequency
          % table at the data level.
 
          Trend   = lower(obj.getOptions('Trend'));
          Sorting = obj.getOptions('SortCategories');
          
          % Step 2: Sort the data according to the criterion in 
          % 'SortCategories'. Ordinal data are treated as numeric data and
          % no sorting is applied
          if ~dbc.isOrdinal
              dbc = obj.applySorting(dbc,Sorting);
          end
          
          FreqTab = dbc.getFrequencyTable('BinLevel',false);          
          
          % Step 4: Compute cutpoint indices, then update the binning
          % container 'dbc'.
          CutPointIndices = obj.runMonotoneAdjacentPoolingMainForCategorical(FreqTab,...
              Trend);
          
          %if ~isempty(CutPointIndices)
              dbc = obj.setBinsByCutPointIndices(dbc,CutPointIndices);
          %end
          
      end
      
      
      function [CutPointIndices, nanBin] = runMonotoneAdjacentPoolingMainForNumeric(~,FreqTab,NanFreqTab,Trend)
          % This is the main MAPA implementation. This function finds the
          % indices where the minima or maxima of the cumulative 
          % characteristic (cumulative Goods) - are located.
          
          Freq    = sum(FreqTab,2);
          NumBins = size(FreqTab,1);
          
          % Step 0: No binning if we have two bins or less.
          % if NumBins <= 2
              % CutPointIndices = [];
              % nanBin = 0;
              % return;
          % end
          
          % Step 1: Set the function associated with the input 'Trend'.
          if strcmpi(Trend,'auto')
              % When no trend is specified by the user, it is set to 'Auto'.
              % This is the default. Given this value, we search for the
              % values iMin and iMax of the index where the minimum and 
              % maximum cumulative Good ratios happen. If iMin <= iMax, 
              % then we expect the trend to be increasing. Otherwise, we 
              % expect it to be decreasing.

              CumGoodRatio = cumsum(FreqTab(:,1)) ./ cumsum(Freq);
              [~,iMax] = max(CumGoodRatio);
              [~,iMin] = min(CumGoodRatio);
              if iMin <= iMax
                  TrendFcn = @min; % Monotone increasing
              else
                  TrendFcn = @min; % Monotone decreasing
              end
          else
              if strcmpi(Trend,'increasing')
                  TrendFcn = @min;
              else
                  TrendFcn = @max;
              end
          end

          % Step 2: Run MAPA (while loop).
          % Note: our convention is that the first column of the frequency 
          % table corresponds to "Goods" and the second column to "Bads".
          CutPointIndices = [];
          iStart  = 1;

          while iStart < NumBins
              CumGood  = cumsum(FreqTab(iStart:end,1));
              CumFreq  = cumsum(Freq(iStart:end));
              CumGoodRatio = CumGood ./ CumFreq;
              Incr  = find(CumGoodRatio==TrendFcn(CumGoodRatio),1,'last');
              iCutPoint = iStart + Incr;
              if iCutPoint <= NumBins
                  CutPointIndices = [CutPointIndices(:)' iCutPoint];
              end
              iStart = iCutPoint; 
          end
          
          totalNoNan = sum(Freq);
          totalNan = sum(NanFreqTab);
          
          totalFreq = totalNoNan + totalNan;
          
          nanGoodRatio = NanFreqTab(1) / totalNan;
          nanBin = 0;
          while true
              if ~isempty(CutPointIndices)
                  firstFreq = sum(Freq(1:(CutPointIndices(1) - 1)));
                  firstGoodCnt = sum(FreqTab(1:(CutPointIndices(1) - 1), 1));
              else
                  firstFreq = sum(Freq);
                  firstGoodCnt = sum(FreqTab(:, 1));
              end
              firstGoodRatio = firstGoodCnt / firstFreq;
              if nanBin == 0 && (totalNan < 0.05 * totalFreq || nanGoodRatio == 0 || nanGoodRatio == 1)
                  if strcmpi(Trend,'increasing')
                      if nanGoodRatio <= firstGoodRatio
                          nanBin = 1;
                      end                  
                  elseif strcmpi(Trend,'decreasing')
                      if nanGoodRatio >= firstGoodRatio
                          nanBin = 1;
                      end
                  end
              end              
              if nanBin == 1
                  firstFreq = firstFreq + totalNan;
                  firstGoodCnt = firstGoodCnt + NanFreqTab(1);
              end
              firstGoodRatio = firstGoodCnt / firstFreq;
              if strcmpi(Trend,'increasing') && (firstFreq < 0.05 * totalFreq || firstGoodRatio == 0)
                      CutPointIndices(1) = [];
              elseif strcmpi(Trend,'decreasing') && (firstFreq < 0.05 * totalFreq || firstGoodRatio == 1)
                      CutPointIndices(1) = [];
              else
                  break;                  
              end                  
          end
          
          while true
              if ~isempty(CutPointIndices)
                  lastFreq = sum(Freq(CutPointIndices(end):end));
                  lastGoodCnt = sum(FreqTab(CutPointIndices(end):end, 1));
              else
                  lastFreq = sum(Freq);
                  lastGoodCnt = sum(FreqTab(:, 1));
              end
              lastGoodRatio = lastGoodCnt / lastFreq;
              if nanBin == 0 && (totalNan < 0.05 * totalFreq || nanGoodRatio == 0 || nanGoodRatio == 1)
                  if strcmpi(Trend,'increasing')
                      if nanGoodRatio >= lastGoodRatio
                          nanBin = 2;
                      end                  
                  elseif strcmpi(Trend,'decreasing')
                      if nanGoodRatio <= lastGoodRatio
                          nanBin = 2;
                      end
                  end
              end
              if nanBin == 2
                  lastFreq = lastFreq + totalNan;
                  lastGoodCnt = lastGoodCnt + NanFreqTab(1);
              end
              lastGoodRatio = lastGoodCnt / lastFreq;
              if strcmpi(Trend,'increasing') && (lastFreq < 0.05 * totalFreq || lastGoodRatio == 1)
                      CutPointIndices(end) = [];
              elseif strcmpi(Trend,'decreasing') && (lastFreq < 0.05 * totalFreq || lastGoodRatio == 0)
                      CutPointIndices(end) = [];
              else
                  break;                  
              end                  
          end
      end
      
      function [CutPointIndices, nanBin] = runMonotoneAdjacentPoolingMainForCategorical(~,FreqTab,Trend)
          % This is the main MAPA implementation. This function finds the
          % indices where the minima or maxima of the cumulative 
          % characteristic (cumulative Goods) - are located.
          
          Freq    = sum(FreqTab,2);
          NumBins = size(FreqTab,1);
          
          % Step 0: No binning if we have two bins or less.
          %if NumBins <= 2
              %CutPointIndices = [];
              %nanBin = 0;
              %return;
          %end
          
          % Step 1: Set the function associated with the input 'Trend'.
          if strcmpi(Trend,'auto')
              % When no trend is specified by the user, it is set to 'Auto'.
              % This is the default. Given this value, we search for the
              % values iMin and iMax of the index where the minimum and 
              % maximum cumulative Good ratios happen. If iMin <= iMax, 
              % then we expect the trend to be increasing. Otherwise, we 
              % expect it to be decreasing.

              CumGoodRatio = cumsum(FreqTab(:,1)) ./ cumsum(Freq);
              [~,iMax] = max(CumGoodRatio);
              [~,iMin] = min(CumGoodRatio);
              if iMin <= iMax
                  TrendFcn = @min; % Monotone increasing
                  Trend = 'increasing';
              else
                  TrendFcn = @min; % Monotone decreasing
                  Trend = 'decreasing';
              end
          else
              if strcmpi(Trend,'increasing')
                  TrendFcn = @min;
              else
                  TrendFcn = @max;
              end
          end

          % Step 2: Run MAPA (while loop).
          % Note: our convention is that the first column of the frequency 
          % table corresponds to "Goods" and the second column to "Bads".
          CutPointIndices = [];
          iStart  = 1;

          while iStart < NumBins
              CumGood  = cumsum(FreqTab(iStart:end,1));
              CumFreq  = cumsum(Freq(iStart:end));
              CumGoodRatio = CumGood ./ CumFreq;
              Incr  = find(CumGoodRatio==TrendFcn(CumGoodRatio),1,'last');
              iCutPoint = iStart + Incr;
              if iCutPoint <= NumBins
                  CutPointIndices = [CutPointIndices(:)' iCutPoint];
              end
              iStart = iCutPoint; 
          end          
          totalFreq = sum(Freq);         
          
          while true
              if ~isempty(CutPointIndices)
                  firstFreq = sum(Freq(1:(CutPointIndices(1) - 1)));
                  firstGoodCnt = sum(FreqTab(1:(CutPointIndices(1) - 1), 1));
              else
                  break;
              end
              firstGoodRatio = firstGoodCnt / firstFreq;              
              if strcmpi(Trend,'increasing') && (firstFreq < 0.05 * totalFreq || firstGoodRatio == 0)
                      CutPointIndices(1) = [];
              elseif strcmpi(Trend,'decreasing') && (firstFreq < 0.05 * totalFreq || firstGoodRatio == 1)
                      CutPointIndices(1) = [];
              else
                  break;                  
              end                  
          end
          
          while true
              if ~isempty(CutPointIndices)
                  lastFreq = sum(Freq(CutPointIndices(end):end));
                  lastGoodCnt = sum(FreqTab(CutPointIndices(end):end, 1));
              else
                  break;
              end
              lastGoodRatio = lastGoodCnt / lastFreq;              
              if strcmpi(Trend,'increasing') && (lastFreq < 0.05 * totalFreq || lastGoodRatio == 1)
                      CutPointIndices(end) = [];
              elseif strcmpi(Trend,'decreasing') && (lastFreq < 0.05 * totalFreq || lastGoodRatio == 0)
                      CutPointIndices(end) = [];
              else
                  break;                  
              end                  
          end
      end
      
      function CutPointIndices = runMonotoneAdjacentPoolingMain(~,FreqTab,Trend)
          % This is the main MAPA implementation. This function finds the
          % indices where the minima or maxima of the cumulative 
          % characteristic (cumulative Goods) - are located.
          
          Freq    = sum(FreqTab,2);
          NumBins = size(FreqTab,1);
          
          % Step 0: No binning if we have two bins or less.
          %if NumBins <= 2
              %CutPointIndices = [];
             % return;
          %end
          
          % Step 1: Set the function associated with the input 'Trend'.
          if strcmpi(Trend,'auto')
              % When no trend is specified by the user, it is set to 'Auto'.
              % This is the default. Given this value, we search for the
              % values iMin and iMax of the index where the minimum and 
              % maximum cumulative Good ratios happen. If iMin <= iMax, 
              % then we expect the trend to be increasing. Otherwise, we 
              % expect it to be decreasing.

              CumGoodRatio = cumsum(FreqTab(:,1)) ./ cumsum(Freq);
              [~,iMax] = max(CumGoodRatio);
              [~,iMin] = min(CumGoodRatio);
              if iMin <= iMax
                  TrendFcn = @min; % Monotone increasing
              else
                  TrendFcn = @min; % Monotone decreasing
              end
          else
              if strcmpi(Trend,'increasing')
                  TrendFcn = @min;
              else
                  TrendFcn = @max;
              end
          end

          % Step 2: Run MAPA (while loop).
          % Note: our convention is that the first column of the frequency 
          % table corresponds to "Goods" and the second column to "Bads".
          CutPointIndices = [];
          iStart  = 1;

          while iStart < NumBins
              CumGood  = cumsum(FreqTab(iStart:end,1));
              CumFreq  = cumsum(Freq(iStart:end));
              CumGoodRatio = CumGood ./ CumFreq;
              Incr  = find(CumGoodRatio==TrendFcn(CumGoodRatio),1,'last');
              iCutPoint = iStart + Incr;
              if iCutPoint <= NumBins
                  CutPointIndices = [CutPointIndices(:)' iCutPoint];
              end
              iStart = iCutPoint; 
          end
          
          
          
      end
      
   end  
      
    
end
