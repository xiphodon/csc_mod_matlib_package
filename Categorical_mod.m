classdef Categorical_mod < Container_mod
% CATEGORICAL Categorical data binning container.
%
%   The categorical container derives from the data binning Container base
%   class, and handles binning use cases for categorical data.
%
%   This class is used for manual and automatic binning applications (see
%   CREDITSCORECARD class for example).
%
%   This is an internal function that is not meant to be called directly
%   by the user.

% Copyright 2014-2015 The MathWorks, Inc. 

   properties (Dependent, SetAccess = private)
      SortedPredictorCategories % Categories observed in data, sorted
      NumBins % Number of bins
   end
   
   properties (Access = protected)
      CatGrouping = table; % Grouping information, from categories to bins
      IsOrdinal = false; % Is the data ordinal
   end
   
   methods
      
      % Constructor
      function obj = Categorical_mod(PredictorData,ResponseData,varargin)
         % CATEGORICAL Create a categorical data binning container.
         % 
         % Inputs are the same as those for the base Container class,
         % except that PredictorData must be a categorical array.
         
         if nargin<2
            ResponseData = [];
         end
         if nargin<1
            PredictorData = [];
         end
         
         if ~iscategorical(PredictorData)
            error(message('finance:internal:finance:binning:Categorical:InvalidPredictorData'));
         end
         
         if size(PredictorData,2)>1 || size(ResponseData,2)>1
            error(message('finance:internal:finance:binning:Categorical:InvalidInputSize'));
         end
         
         obj = obj@Container_mod(PredictorData,ResponseData,varargin{:});
         
         obj.CatGrouping.Category = obj.PredictorLabels;
         obj.CatGrouping.BinNumber = (1:length(obj.PredictorLabels))';
         obj.IsOrdinal = isordinal(obj.PredictorData);
         
      end
      
      % Public properties
      function cats = get.SortedPredictorCategories(obj)
         cats = obj.CatGrouping.Category;
      end
      
      % Dependent properties
      function n = get.NumBins(obj)
         ug = unique(obj.CatGrouping.BinNumber);
         n = length(ug);
         if ismember(0,ug)
            n = n-1;
         end
      end
      
      function values = getCatGrouping(obj)
         values = obj.CatGrouping;
      end
      
      function obj = setCatGrouping(obj,CatGrouping)
         % SETCATGROUPING Modify bins according to the given grouping.
         %
         % The binning rules for categorical are defined via grouping. The
         % CatGrouping, or category grouping, input is a table with two
         % columns: Category and BinNumber. Each category in column one
         % falls under the bin number in column two.
         %
         % Notes:
         %   o User-defined bin labels are overridden each time a new
         %     category grouping is provided. This is because bin labels
         %     often contain references to the binning rules, and changing
         %     binning rules may cause labels and binning rules to be out
         %     of sync.
         
         if ~istable(CatGrouping)
            error(message('finance:internal:finance:binning:Categorical:InvalidCatGrouping'));
         end
         
         [m,n] = size(CatGrouping);
         
         CGVarNames = CatGrouping.Properties.VariableNames;
         if n ~= 2 || ~all(ismember({'Category','BinNumber'},CGVarNames))
            error(message('finance:internal:finance:binning:Categorical:InvalidCatGroupingCols'));
         end
         
         CatsObj = obj.CatGrouping.Category;
         CatsInput = CatGrouping.Category;
         if m < obj.NumRows || ~all(ismember(CatsObj,CatsInput))
            error(message('finance:internal:finance:binning:Categorical:InvalidCatGroupingCats'));
         end
         
         Grouping = CatGrouping.BinNumber;
         if ~isnumeric(Grouping) || ~all(mod(Grouping,1)==0)
            error(message('finance:internal:finance:binning:Categorical:InvalidCatGroupingBins'));
         end
         
         if any(Grouping<0)
            error(message('finance:internal:finance:binning:Categorical:NegativeCatGroupingBins'));
         end
         
         if obj.IsOrdinal
            if any(diff(Grouping)<0)
               warning(message('finance:internal:finance:binning:Categorical:NonMonotonicGrouping'));
               obj.IsOrdinal = false;
            end
         end
         
         % Make sure bin numbers are consecutive
         uGrouping = unique(Grouping);
         uGrouping(uGrouping == 0) = []; % do not reset out of range
         NumGroups = length(uGrouping);
         
         for ii = 1:NumGroups
            ind = (Grouping == uGrouping(ii));
            CatGrouping.BinNumber(ind) = ii;
         end

         obj.CatGrouping = CatGrouping;
         
         % Reset user labels to empty, since they may now be out of sync
         obj.UserLabels = {};
         
      end
      
      function flag = isOrdinal(obj)
         flag = obj.IsOrdinal;
      end
      
      function labels = getLabels(obj)
         % GETLABELS Get bin labels.
         %
         % For categorical data, bin labels are simply 'Group1', 'Group2',
         % etc.

         if ~isempty(obj.UserLabels)
            labels = obj.UserLabels;
            return;
         end
         
         if obj.NumBins <= 0
            labels = {};
            return;
         end
         
         CatLabels = obj.CatGrouping.Category;
         Grouping = obj.CatGrouping.BinNumber;
         labels = cell(obj.NumBins,1);
                 
         % Remove out of range
         IndOutOfRange = (Grouping == 0);
         Grouping(IndOutOfRange)  = [];
         CatLabels(IndOutOfRange) = [];
         
         if obj.NumBins == length(CatLabels)
            [~,InvGrouping] = sort(Grouping);
            labels = CatLabels(InvGrouping);
            return;
         end
         
         for ii=1:obj.NumBins
            labels{ii} = ['Group' num2str(ii)];
         end
         
      end
      
      function members = getMembers(obj,varargin)
         % GETMEMBERS Get bin members.
         %
         % This method returns the categories in a given bin. The output is
         % a cell array, with one row per bin. The element in each row can
         % be in turn a cell array of strings with the category names, or
         % it can be a single string where all the category names in a bin
         % are concatenated.
         
         p = inputParser;
         p.addParameter('Format','Cell',@(x)ismember(lower(x),{'cell','string'}));
         p.parse(varargin{:});
         Format = p.Results.Format;
         
         if obj.NumBins <= 0
            members = {};
            return;
         end
         
         members = cell(obj.NumBins,1);
         ugrouping = unique(obj.CatGrouping.BinNumber);
         ugrouping(ugrouping==0) = []; % drop out of range
         for ii = 1:obj.NumBins
            ind = obj.CatGrouping.BinNumber==ugrouping(ii);
            members{ii} = obj.CatGrouping.Category(ind)';
         end

         if strcmpi(Format,'string')
            for ii=1:obj.NumBins
               s = sprintf('%s,',members{ii}{:});
               s(end) = '';
               members{ii} = s;
            end
         end
         
      end
      
      function [obj,ind,indinv] = sortPredictorCategories(obj,varargin)
         % SORTPREDICTORCATEGORIES Sort the list of categories in the data.
         %
         % The categorical data has an associated list of possible
         % categories, and the order in this list usually is not determined
         % by the user. However, this ordering is used to, for example,
         % display the frequency table. For some binning use cases,
         % re-ordering may make a difference in the final outcome of a
         % procedure or algorithm. This method allows clients to re-order
         % the categories.
         %
         % It supports a few sorting criteria such as ordering by row (or
         % histogram) count, by odds, or by frequencies in a particular
         % column. It allows sorting in ascending or descending mode. It
         % also allows the caller to pass a vector of sorted indices.
         % 
         
         p = inputParser;
         
         p.addParameter('SortedIndices',[],@(x)all((isnumeric(x)&(x>0)&(mod(x,1)==0))));
         p.addParameter('Criterion',{'RowCount'},...
            @(x)all(ismember(lower(x),...
            {'rowcount','odds','colcount'})));
         p.addParameter('Mode',{'Ascending'},...
            @(x)all(ismember(lower(x),...
            {'ascending','descending'})));
         p.addParameter('ColNumber',[],...
            @(x)isnumeric(x)&(all(isnan(x)|(x>=1&x<=obj.NumCols&mod(x,1)==0))));
         
         p.parse(varargin{:});

         ind = p.Results.SortedIndices;
         if ~isempty(ind)
            [obj,indinv] = obj.sortPredictorCategoriesByIndices(ind);
            return;
         end
         
         Criterion = p.Results.Criterion;
         if ischar(Criterion)
            Criterion = {Criterion};
         end
         
         Mode = p.Results.Mode;
         if ischar(Mode)
            Mode = {Mode};
         end
         
         if length(Mode)==1
            Mode = repmat(Mode,size(Criterion));
         elseif length(Mode)>1
            if length(Mode)~=length(Criterion)
               error(message('finance:internal:finance:binning:Categorical:InvalidSortingMode'));
            end
         end
         
         ColNumber = p.Results.ColNumber;
         
         ColCountInd = strcmpi('colcount',Criterion);
         NumColCount = sum(ColCountInd);
         if isempty(ColNumber)
            if NumColCount == 1
               ColNumber = NaN(size(Criterion));
               ColNumber(ColCountInd) = obj.NumCols;
            elseif NumColCount > 1
               error(message('finance:internal:finance:binning:Categorical:MissingClassNumber'));
            end % else if NumCountClass == 0, nothing to do
         else
            if length(ColNumber)~=length(Criterion)
               error(message('finance:internal:finance:binning:Categorical:InvalidColNumber'));
            else
               NotNaNInd = ~isnan(ColNumber);
               if any(NotNaNInd~=ColCountInd)
                  error(message('finance:internal:finance:binning:Categorical:InvalidColCountInd'));
               end
            end
         end

         if ismember('Odds',Criterion)&&obj.NumCols~=2
            error(message('finance:internal:finance:binning:Categorical:MismatchOddsNumCols'));
         end
         
         NumCrit = length(Criterion);
         FreqTab = obj.FrequencyTable;
         SortMatrix = zeros(obj.NumRows,NumCrit);
         SortFlag = 1:NumCrit;
         
         for ii = 1:NumCrit
            switch lower(Criterion{ii})
               case 'rowcount'
                  SortMatrix(:,ii) = sum(FreqTab,2);
               case 'odds'
                  SortMatrix(:,ii) = FreqTab(:,1)./FreqTab(:,2);
               case 'colcount'
                  SortMatrix(:,ii) = FreqTab(:,ColNumber(ii));
            end
            if strcmpi(Mode{ii},'descending')
               SortFlag(ii) = -SortFlag(ii);
            end
         end
         
         [~,ind] = sortrows(SortMatrix,SortFlag);

         [obj,indinv] = obj.sortPredictorCategoriesByIndices(ind);

      end
      
      function bdata = binData(obj,varargin)
         % BINDATA Bin given data into predetermined bins.
         %
         % Given the binning rules given in the category grouping, this
         % method can map observed data into bins. By default it returns a
         % numeric vector with bin numbers. The user can pass either labels
         % or numeric values corresponding to each bin, and this method can
         % then map the data into those values, instead of bin numbers.
         % One extra label or value can be used for missing or out-of-range
         % data.
         
         p = inputParser;
         p.addOptional('data',obj.PredictorData,@iscategorical);
         p.addParameter('Values',[],@(x)isnumeric(x)||(iscell(x)&&all(cellfun(@ischar,x))));
         p.parse(varargin{:});
         
         data = p.Results.data(:);
         Values = p.Results.Values(:);
         
         if isempty(data)
            if obj.IsCompact
               error(message('finance:internal:finance:binning:Categorical:EmptyBinDataCompact'));
            else
               error(message('finance:internal:finance:binning:Categorical:EmptyBinData'));
            end
         end
         
         NumValues = length(Values);
         if ~isempty(Values)&&~(NumValues==obj.NumBins || NumValues==obj.NumBins+1)
            error(message('finance:internal:finance:binning:Categorical:InvalidBinDataValues'));
         end

         if ~isempty(Values)&&iscell(Values)
            labels = Values(1:obj.NumBins);
         else
            labels = obj.getLabels;
         end
         isord = obj.IsOrdinal;
         cg = obj.getCatGrouping;
         bdata = categorical(data,cg.Category,labels(cg.BinNumber),'Ordinal',isord);
         bdata = reordercats(bdata,labels);
         
         if isempty(Values) % return bin number
            % Out of range are automatically converted to NaN
            bdata = double(bdata);
         else
            if isnumeric(Values) % Return numeric values
               bdata = double(bdata);
               nanind = isnan(bdata);
               bdata(~nanind) = Values(bdata(~nanind));
               if NumValues == obj.NumBins+1
                  bdata(nanind) = Values(end);
               end
            else % handling of missing value, categorical
               if NumValues == obj.NumBins+1
                  nanind = isnan(double(bdata));
                  bdata = addcats(bdata,Values{end},'After',Values{end-1});
                  bdata(nanind) = Values{end};
               end
            end         
         end
         
      end
       
      function Stats = getPredictorStatistics(obj)
          % Returns a table, 'Stats', with the row names being the names of 
          % the categories, with corresponding counts.

          Categories = obj.SortedPredictorCategories;
          Count = sum(obj.getFrequencyTable('BinLevel',false),2);
              
          Stats = table(Count,'RowNames',Categories);

      end
      
      
   end

   methods (Access = protected)
      
      function [FreqTab, NanFreqTab] = getBinLevelFrequencyTable(obj)
         % GETBINLEVELFREQUENCYTABLE Return frequencies at bin level.
         %
         % Compute the frequency table at the bin level.
         NanFreqTab = obj.nanFrequencyTable;
         if obj.NumBins <= 0
            FreqTab = [];
         else
            NumClasses = size(obj.FrequencyTable,2);
            FreqTab = zeros(obj.NumBins,NumClasses);
            ugrouping = unique(obj.CatGrouping.BinNumber);
            ugrouping(ugrouping==0) = []; % drop out of range
            for ii = 1:obj.NumBins
               ind = obj.CatGrouping.BinNumber==ugrouping(ii);
               FreqTab(ii,:) = sum(obj.FrequencyTable(ind,:),1);
            end
         end

      end
      
      function [obj,SortIndInv] = sortPredictorCategoriesByIndices(obj,SortedIndices)
         % SORTPREDICTORCATEGORIESBYINDICES Sort data categories according
         %   to a given re-ordering vector of indices.
         %
         % This is a companion to sortPredictorCategories. This method
         % re-orders the categories based on a vector of indices.
         % sortPredictorCategories supports re-ordering based on multiple
         % criteria, whereas this method can only take a vector of sorted
         % indices.
         
         cg = obj.CatGrouping;
         NumCats = size(cg,1);
         InRangeInd = obj.CatGrouping.BinNumber>0;
         OutOfRangeInd = ~InRangeInd;
         NumInRange = sum(InRangeInd);
         
         if length(SortedIndices)~=NumInRange
            error(message('finance:internal:finance:binning:Categorical:InvalidSortedIndices'));
         end
         
         FreqTabInRange = obj.FrequencyTable(InRangeInd,:);         
         FreqTabOutOfRange = obj.FrequencyTable(OutOfRangeInd,:);
         obj.FrequencyTable = [FreqTabInRange;FreqTabOutOfRange];
         obj.FrequencyTable(1:NumInRange,:) = FreqTabInRange(SortedIndices,:);

         LabelsInRange = obj.PredictorLabels(InRangeInd);
         LabelsOutOfRange = obj.PredictorLabels(OutOfRangeInd);
         obj.PredictorLabels = [LabelsInRange;LabelsOutOfRange];
         obj.PredictorLabels(1:NumInRange) = LabelsInRange(SortedIndices);
         
         OriginalOrder = all(cg.BinNumber==(1:NumCats)');
         
         cgInRange = cg(InRangeInd,:);
         cgOutOfRange = cg(OutOfRangeInd,:);
         cg = [cgInRange;cgOutOfRange];
         cg(1:NumInRange,:) = cgInRange(SortedIndices,:);
         
         if OriginalOrder
            cg.BinNumber = (1:NumCats)';
         end
         
         obj = obj.setCatGrouping(cg);
         obj.IsOrdinal = false;
         
         [~,SortIndInv] = sort(SortedIndices);
         
      end
      
   end
   
end
