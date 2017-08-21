classdef Numeric_mod < Container_mod
% NUMERIC Numeric data binning container.
%
%   The numeric container derives from the data binning Container base
%   class, and handles binning use cases for numeric data.
%
%   This class is used for manual and automatic binning applications (see
%   CREDITSCORECARD class for example).
%
%   This is an internal function that is not meant to be called directly
%   by the user.

% Copyright 2014-2015 The MathWorks, Inc. 

   properties (SetAccess = private)
      MinValue = -Inf; % Minimum value (not necessarily observed in the data)
      MaxValue = Inf;  % Minimum value (not necessarily observed in the data)
      SortedPredictorValues = []; % Values observed in data, sorted
   end
   
   properties (Dependent, SetAccess = private)
      NumBins % Number of bins
   end
   
   properties (Dependent, Access = protected)
      NumBinsInRange % Bin with data in range (>=MinVal, <=MaxVal)
      FirstPointer   % Pointer to the beginning of the first bin
      LastPointer    % Pointer to the end of last bin, plus one
      HasOutOfRange  % If true, there are observations are out of range
   end
   
   properties (Access = protected)
      BinPointers = [];  % Pointers to the beginning of bins (one extra at the end)
      CutPoints = [];    % Cut points determining the bins
      IsInteger = false; % Is the data integer
      BasicStatistics    % Table with basic predictor statistics (min, max, etc.)
      NanBin = 0;
      HasNaN = false;
   end
   
   methods
      
      % Constructor
      function obj = Numeric_mod(PredictorData,ResponseData,varargin)
         % NUMERIC Create a numeric data binning container.
         % 
         % Inputs are the same as those for the base Container class,
         % except that PredictorData must be a numeric array.
         
         % Set the warning backtrace option to 'off' to avoid the command
         % line stack pointing to the methods of the internal package.         
         %warning off backtrace
         if nargin<2
            ResponseData = [];
         end
         if nargin<1
            PredictorData = [];
         end
         
         if ~isnumeric(PredictorData)
            error(message('finance:internal:finance:binning:Numeric:InvalidPredictorData'));
         end
         
         if size(PredictorData,2)>1 || size(ResponseData,2)>1
            error(message('finance:internal:finance:binning:Numeric:InvalidInputSize'));
         end
         
         obj = obj@Container_mod(PredictorData,ResponseData,varargin{:});
         
         obj.HasNaN = any(isnan(PredictorData));    
         obj.SortedPredictorValues = cellfun(@str2num,obj.PredictorLabels);

         obj.BinPointers = obj.FirstPointer:obj.LastPointer;
         obj.CutPoints = obj.SortedPredictorValues(obj.BinPointers(2:end-1));
         obj.IsInteger = all(mod(obj.SortedPredictorValues,1)==0);
         
         % Get basic predictor statistics. This is called here so that
         % these statistics are available even if the container is compact.
         % However, because they are computed at construction, changes to
         % MinValue and MaxValue will not affect the basic statistics. This
         % is just a convenience for an initial summary of the data, the
         % binning operations are unrelated to these basic statistics.
         obj = getBasicStatistics(obj,PredictorData);
         
      end
      
      % Dependent properties
      function n = get.NumBins(obj)
         n = length(obj.BinPointers)-1;
         n = max(n,0);
      end
      
      function n = get.FirstPointer(obj)
         n = find(obj.SortedPredictorValues>=obj.MinValue,1,'first');
         if isempty(n)
            n = obj.NumRows+1;
         end
      end

      function n = get.LastPointer(obj)
         n = find(obj.SortedPredictorValues<=obj.MaxValue,1,'last')+1;
         if isempty(n)
            n = 1;
         end
      end

      function flag = get.HasOutOfRange(obj)
         flag = any(obj.SortedPredictorValues<obj.MinValue |...
            obj.SortedPredictorValues>obj.MaxValue);
      end

      function values = getCutPoints(obj)
         values = obj.CutPoints;
      end
      
      function obj = setCutPoints(obj,CutPoints, NanBin)
         % SETCUTPOINTS Modify bins according to the given cut points.
         %
         % The binning rules for numeric are defined via cut points. The
         % cut points C1, ..., CN define N+1 bins, with:
         %
         %      Bin 1: [MinValue,C1)
         %      Bin 2: [C1,C2)
         %      ...
         %      Bin N+1: [CN,MaxValue]
         %
         % Notes:
         %   o User-defined bin labels are overridden each time a new
         %     category grouping is provided. This is because bin labels
         %     often contain references to the binning rules, and changing
         %     binning rules may cause labels and binning rules to be out
         %     of sync.
         obj.NanBin = NanBin;
         if isempty(CutPoints)
            obj.BinPointers = [obj.FirstPointer;obj.LastPointer];
            obj.CutPoints = CutPoints;
            return;
         end
         
         if ~isnumeric(CutPoints) || any(diff(CutPoints)<0)
            error(message('finance:internal:finance:binning:Numeric:InvalidCutPoints'));
         end
         
         if any(diff(CutPoints)==0)
            warning(message('finance:internal:finance:binning:Numeric:DuplicateCutPoints'));
            CutPoints = unique(CutPoints);
         end
         
         if obj.MinValue >= min(CutPoints)
            warning(message('finance:internal:finance:binning:Numeric:CutPointsBelowMin'));
            CutPoints = CutPoints(obj.MinValue < CutPoints);
         end
         
         if obj.MaxValue < max(CutPoints)
            warning(message('finance:internal:finance:binning:Numeric:CutPointsAboveMax'));
            CutPoints = CutPoints(obj.MaxValue >= CutPoints);
         end
         
         if any(isnan([obj.MinValue;obj.MaxValue;CutPoints(:)]))
             error(message('finance:internal:finance:binning:Numeric:NaNCutPoints'));
         end
         
         obj.CutPoints = CutPoints(:);
         
         SortedValues = obj.SortedPredictorValues;
         obj.BinPointers = zeros(length(CutPoints)+2,1);
         
         BelowRange = CutPoints(:)<min(obj.SortedPredictorValues);
         NumBelow = sum(BelowRange);
         obj.BinPointers(1:1+NumBelow) = obj.FirstPointer;
         
         AboveRange = CutPoints(:)>max(obj.SortedPredictorValues);
         NumAbove = sum(AboveRange);
         obj.BinPointers(end-NumAbove:end) = obj.LastPointer;
         
         InRange = true(size(CutPoints(:)))&~BelowRange&~AboveRange;
         NumIn = sum(InRange);
         InRangeCutPoints = CutPoints(InRange);
         for ii=1:NumIn
            obj.BinPointers(1+NumBelow+ii)=...
               find(SortedValues>=InRangeCutPoints(ii),1,'first');
         end
         
         % Reset user labels to empty, since they may now be out of sync
         obj.UserLabels = {};
         
      end

      function labels = getLabels(obj)
         % GETLABELS Get bin labels.
         %
         % For numeric data, bin labels are given in terms of the cut
         % points. If the cut points are C1, ..., CN, the labels for the
         % corresponding N+1 bins are
         %
         %      Bin 1: '[MinValue,C1)'
         %      Bin 2: '[C1,C2)'
         %      ...
         %      Bin N+1: '[CN,MaxValue]'
         %
         
         if ~isempty(obj.UserLabels)
            labels = obj.UserLabels;
            return;
         end
         
         if obj.NumBins <= 0
            labels = {};
            return;
         end
         
         if obj.NumBins == length(obj.SortedPredictorValues) && ~obj.HasOutOfRange
            labels = obj.PredictorLabels;
            if obj.HasNaN
                if obj.NanBin == 1
                    labels{1} = [labels{1} ' & Missing Value'];
                elseif obj.NanBin == 2
                    labels{end} = [labels{end} ' & Missing Value'];
                else
                    labels{end + 1} = 'Missing Value';   
                end
            end
            return;
         end
         
         cp = obj.getCutPoints;

         labels = cell(obj.NumBins,1);
         
         if obj.NumBins == 1
            labels{1} = ['[' num2str(obj.MinValue) ',' num2str(obj.MaxValue) ']'];
         else
            % First label
            labels{1} = ['[' num2str(obj.MinValue) ',' num2str(cp(1)) ')'];
            % Middle labels
            for ii=2:obj.NumBins-1
               labels{ii} = ['[' num2str(cp(ii-1)) ',' num2str(cp(ii)) ')'];
            end
            % Last label
            labels{obj.NumBins} = ['[' num2str(cp(obj.NumBins-1)) ',' num2str(obj.MaxValue) ']'];
         end
         if obj.HasNaN
             if obj.NanBin == 1
                labels{1} = [labels{1} ' & Missing Value'];
             elseif obj.NanBin == 2
                labels{end} = [labels{end} ' & Missing Value'];
             else
                labels{end + 1} = 'Missing Value';   
             end  
         end
      end

      function obj = setMinValue(obj,minvalue)
         % SETMINVALUE Set minimum value for data.
         %
         % The range of values observed in the sample may not include the
         % true minimum or maximum possible values for this type of data.
         % This method is used to set a minimum value for the data.
         
         if ~isnumeric(minvalue) || ~isscalar(minvalue)
            error(message('finance:internal:finance:binning:Numeric:InvalidMinValue'));
         end
         
         if isnan(minvalue)
             error(message('finance:internal:finance:binning:Numeric:NaNCutPoints'));
         end
         
         if minvalue > obj.MaxValue
            error(message('finance:internal:finance:binning:Numeric:MinAboveMax'));
         end
         
         if ~isempty(obj.SortedPredictorValues) &&...
               minvalue > min(obj.SortedPredictorValues)
            warning(message('finance:internal:finance:binning:Numeric:MinAboveData'));
         end
         
         obj.MinValue = minvalue;

         % Refresh cut points and bin pointers
         cp = obj.getCutPoints;
         if ~isempty(cp) && minvalue >= min(cp)
            warning(message('finance:internal:finance:binning:Numeric:MinAboveCutPoints'));
            cp = cp(minvalue < obj.getCutPoints);
         end
         obj = obj.setCutPoints(cp);
         
      end
      
      function obj = setMaxValue(obj,maxvalue)
         % SETMAXVALUE Set maximum value for data.
         %
         % The range of values observed in the sample may not include the
         % true minimum or maximum possible values for this type of data.
         % This method is used to set a maximum value for the data.
         
         if ~isnumeric(maxvalue) || ~isscalar(maxvalue)
            error(message('finance:internal:finance:binning:Numeric:InvalidMaxValue'));
         end
         
         if isnan(maxvalue)
             error(message('finance:internal:finance:binning:Numeric:NaNCutPoints'));
         end
         
         if maxvalue < obj.MinValue
            error(message('finance:internal:finance:binning:Numeric:MaxBelowMin'));
         end
         
         if ~isempty(obj.SortedPredictorValues) &&...
               maxvalue < max(obj.SortedPredictorValues)
            warning(message('finance:internal:finance:binning:Numeric:MaxBelowData'));
         end
         
         obj.MaxValue = maxvalue;

         % Refresh cut points and bin pointers
         cp = obj.getCutPoints;
         if ~isempty(cp) && maxvalue < max(cp)
            warning(message('finance:internal:finance:binning:Numeric:MaxBelowCutPoints'));
            cp = cp(maxvalue >= obj.getCutPoints);
         end
         obj = obj.setCutPoints(cp);
         
      end
      
      function bdata = binData(obj,varargin)
         % BINDATA Bin given data into predetermined bins.
         %
         % Given the binning rules given in the cut points, this method can
         % map observed data into bins. By default it returns a numeric
         % vector with bin numbers. The user can pass either labels or
         % numeric values corresponding to each bin, and this method can
         % then map the data into those values, instead of bin numbers. One
         % extra label or value can be used for missing or out-of-range
         % data.
         
         p = inputParser;
         p.addOptional('data',obj.PredictorData,@isnumeric);
         p.addParameter('Values',[],@(x)isnumeric(x)||(iscell(x)&&all(cellfun(@ischar,x))));
         p.parse(varargin{:});
         
         data = p.Results.data(:);
         Values = p.Results.Values(:);
         
         if isempty(data)
            if obj.IsCompact
               error(message('finance:internal:finance:binning:Numeric:EmptyBinDataCompact'));
            else
               error(message('finance:internal:finance:binning:Numeric:EmptyBinData'));
            end
         end
         
         NumValues = length(Values);
         if ~isempty(Values)&&~(NumValues==obj.NumBins || NumValues==obj.NumBins+1 || NumValues==obj.NumBins+2)
            error(message('finance:internal:finance:binning:Numeric:InvalidBinDataValues'))
         end

         opts = {};
         if ~isempty(Values)
            if isnumeric(Values)
               opts = {Values(1:obj.NumBins)};
            else % cell, categorical output
               opts = {'categorical',Values(1:obj.NumBins)};
            end
         end
         
         BinRanges = [obj.MinValue; obj.getCutPoints; obj.MaxValue];
         
         bdata = discretize(data,BinRanges,opts{:});
         if obj.HasNaN
             if obj.NanBin == 0
                 bdata(isnan(bdata)) = Values(obj.NumBins + 1);
             elseif obj.NanBin == 1
                 bdata(isnan(bdata)) = Values(1);
             elseif obj.NanBin == 2
                 bdata(isnan(bdata)) = Values(obj.NumBins);
             end
         else 
             if NumValues == obj.NumBins+1
                if isnumeric(Values)
                   bdata(isnan(bdata)) = Values(end);
                else % cell, categorical output
                   bdata = addcats(bdata,Values{end},'After',Values{end-1});
                   bdata(isnan(double(bdata))) = Values{end};
                end
             end
         end
         
      end
      
      function Stats = getPredictorStatistics(obj)
          % Returns a table, 'Stats', with the following row names:
          %
          %   - Min : Minimum value in the predictor data.
          %
          %   - Max : Maximum value in the predictor data.
          %
          %   - Mean: Mean value in the sample.
          %
          %   - Std : Standard deviation of the sample.
          %
          %   The corresponding values are in the 'Values' column.
          %
          
          Stats = obj.BasicStatistics;
          
      end
      
   end
   
   methods (Access = protected)
      
      function [FreqTab, NanFreqTab] = getBinLevelFrequencyTable(obj)
         % GETBINLEVELFREQUENCYTABLE Return frequencies at bin level.
         %
         % Compute the frequency table at the bin level.
         
         if obj.NumBins <= 0
            FreqTab = [];
            NanFreqTab = [];
         else
            NumClasses = size(obj.FrequencyTable,2);
            FreqTab = zeros(obj.NumBins,NumClasses);
            for ii = 1:obj.NumBins
               BinStart = obj.BinPointers(ii);
               BinEnd = obj.BinPointers(ii+1)-1;
               FreqTab(ii,:) = sum(obj.FrequencyTable(BinStart:BinEnd,:),1);
            end
            NanFreqTab = [];
            if obj.NanBin == 1
                FreqTab(1,:) = FreqTab(1,:) + obj.nanFrequencyTable;
            elseif obj.NanBin == 2
                FreqTab(end,:) = FreqTab(end,:) + obj.nanFrequencyTable;
            else
                NanFreqTab = obj.nanFrequencyTable;
            end
         end
      end
      
      function obj = getBasicStatistics(obj,data)
         % GETBASICSTATISTICS Compute basic statistics of predictor data.
         %
         % Compute basic statistics of the predictor, and returns them in
         % the BasicStatistics property of the object. This property is a
         % table. The row names correspond to the statistics, and the
         % values are stored in the 'Value' column.
         %
         %   Min: Minimum value.
         %   Max: Maximum value.
         %   Mean: Mean of the data.
         %   Std: Standard deviation of the data.
         %
         % Note:
         %   For data types other than double or single, numeric precision
         %   may be lost for the standard deviation. The underlying
         %   function 'std' only accepts double or single types. All other
         %   data types are cast as 'double' before computing the standard
         %   deviation.
         
          Min = min(data);
          Max = max(data);
          Mean = mean(data);
          if strcmpi(class(data),'double') || strcmpi(class(data),'single')
             Std = std(data);
          else
             Std = std(double(data));
          end
          
          Value = [Min;Max;Mean;Std];
          RowNames  = {'Min','Max','Mean','Std'};
          obj.BasicStatistics = table(Value,'RowNames',RowNames);

      end
      
   end
   
end
