classdef Container_mod
% CONTAINER Base class for data binning containers.
%
%   The data binning containers are data wrappers specialized for data
%   binning use cases. Their most important property is the frequency table
%   of the data. Their methods support frequency reports, modification of
%   bins, re-labeling of bins, etc. These containers can also work with
%   automatic binning algorithms in the binning package.
%
%   These classes are used for manual and automatic binning applications
%   (see CREDITSCORECARD class for example).
%
%   This is an internal function that is not meant to be called directly
%   by the user.

% Copyright 2014 The MathWorks, Inc. 
   
   properties (SetAccess = private)
      PredictorData = []; % Copy of the predictor data
      PredictorVar = '';  % Predictor's name
      ResponseData = [];  % Copy of the response data
      ResponseVar = '';   % Response name
   end

   properties (Dependent, SetAccess = private)
      HasResponseData % Hase response data (suitable for supervised algorithms)
      IsCompact       % If true, stores frequency table only, not original data
   end
   
   properties (Access = protected)
      FrequencyTable = [];  % Frequency table of the data
      nanFrequencyTable = [];
      NumCols = 0;          % Num cols frequency table
      NumRows = 0;          % Num rows frequency table
      PredictorLabels = {}; % Sorted predictor values or categories in string format
      ResponseLabels = {};  % Response values or categories in string format
      UserLabels = {};      % User-provided bin labels
   end
   
   methods
      
      % Constructor
      function obj = Container_mod(PredictorData,ResponseData,varargin)
         % CONTAINER Create a data binning container.
         % 
         % Inputs:
         %
         %   o PredictorData: A numeric or categorical array.
         %
         %   o ResponseData: A numeric or categorical array of the same
         %      size as PredictorData containing the response (or class)
         %      information corresponding to the predictor data.
         %
         % Optional Inputs:
         %
         %   'Compact'       Flag to indicate whether to keep a copy of the
         %                   original data, or only the frequency table. By
         %                   default it is false, and PredictorData and
         %                   ResponseData keep copies of the original data.
         %
         %   'ResponseOrder' Numeric or cell array containing a sorted list
         %                   of all the response values. The frequency
         %                   table's columns are sorted accordingly. This
         %                   is particularly important for binary responses
         %                   when statistics on the frequencies (e.g.,
         %                   odds, weight of evidence) depend on the order
         %                   of the columns of the frequency table.
         
         p = inputParser;
         p.addParameter('PredictorVar','',@ischar);
         p.addParameter('ResponseVar','',@ischar);
         p.addParameter('Compact',false,@islogical);
         p.addParameter('ResponseOrder',{},@(x)(isnumeric(x)||iscell(x)));
         
         parse(p,varargin{:});
         
         obj.PredictorVar = p.Results.PredictorVar;
         obj.ResponseVar = p.Results.ResponseVar;
         Compact = p.Results.Compact;
         ResponseOrder = p.Results.ResponseOrder;
         
         if ~isnumeric(PredictorData) && ~iscategorical(PredictorData)
            error(message('finance:internal:finance:binning:Container:InvalidPredictorData'));
         end
         
         if isempty(PredictorData)
            error(message('finance:internal:finance:binning:Container:EmptyPredictorData'));
         end
         
         if ~isnumeric(ResponseData) && ~iscategorical(ResponseData)
            error(message('finance:internal:finance:binning:Container:InvalidResponseData'));
         end
         
         if ~Compact
            obj.PredictorData = PredictorData;
            obj.ResponseData = ResponseData;
         end

         if ~isempty(ResponseData)
            [freq,~,~,labels] = crosstab(PredictorData,ResponseData);
            if isnumeric(PredictorData)                
                nanfreq = [sum(isnan(PredictorData) & ResponseData == 0) sum(isnan(PredictorData) & ResponseData == 1)];
            else
                nanfreq = [];
            end
         else
            [freq,~,~,labels] = crosstab(PredictorData);
         end

         obj.FrequencyTable = freq;
         obj.nanFrequencyTable = nanfreq;
         [obj.NumRows,obj.NumCols] = size(freq);
         obj.PredictorLabels = labels(:,1);
         obj.PredictorLabels(cellfun(@isempty,obj.PredictorLabels))=[];
         if ~isempty(ResponseData)
            obj.ResponseLabels = labels(:,2);
            obj.ResponseLabels(cellfun(@isempty,obj.ResponseLabels))=[];
            if ~isempty(ResponseOrder)
               if iscell(ResponseOrder) && any(cellfun(@isnumeric,ResponseOrder))
                  error(message('finance:internal:finance:binning:Container:InvalidResponseOrder'));
               end
               if (isnumeric(ResponseOrder) && iscategorical(ResponseData))
                  error(message('finance:internal:finance:binning:Container:ResponseOrderMismatch1'));
               end
               if (iscell(ResponseOrder) && isnumeric(ResponseData))
                  error(message('finance:internal:finance:binning:Container:ResponseOrderMismatch2'));
               end
               if isnumeric(ResponseOrder)
                  TableOrder = cellfun(@str2num,obj.ResponseLabels);
               else
                  TableOrder = lower(obj.ResponseLabels);
                  ResponseOrder = lower(ResponseOrder(:));
               end
               [ok,ind] = ismember(TableOrder,ResponseOrder);
               if ~all(ok)
                  error(message('finance:internal:finance:binning:Container:ResponseOrderMismatch3'));
               end
               obj.FrequencyTable = obj.FrequencyTable(:,ind);
               if ~isempty(obj.nanFrequencyTable)
                    obj.nanFrequencyTable = obj.nanFrequencyTable(ind);
               end
               obj.ResponseLabels = obj.ResponseLabels(ind);
            end
         end

      end

      % Public properties
      function data = get.PredictorData(obj)
         data = obj.PredictorData;
      end
      
      function labels = get.PredictorLabels(obj)
         labels = obj.PredictorLabels;
      end
      
      function name = get.PredictorVar(obj)
         name = obj.PredictorVar;
      end
      
      function data = get.ResponseData(obj)
         data = obj.ResponseData;
      end

      function labels = get.ResponseLabels(obj)
         labels = obj.ResponseLabels;
      end
      
      function name = get.ResponseVar(obj)
         name = obj.ResponseVar;
      end
      
      % Dependent properties
      function flag = get.HasResponseData(obj)
         flag = obj.NumCols>1;
      end
      
      function flag = get.IsCompact(obj)
         flag = isempty(obj.PredictorData)&&~isempty(obj.FrequencyTable);
      end
      
      % General methods
      function [FreqTab, NanFreqTab] = getFrequencyTable(obj,varargin)
         % GETFREQUENCYTABLE Return frequency table of given data.
         %
         % By default, getFrequencyTable returns the frequency table of the
         % binned data (if data has already been binned). However, it
         % supports additional outputs.
         %
         % Optional Inputs:
         %
         % 'BinLevel'  Flag to indicate whether the frequencies of interest
         %             are at bin level (true, default value), or at the
         %             observation level (if false).
         %
         % 'Normalize' Use this option to indicate whether the frequency
         %             table should be normalized. Options are:
         %
         %             o 'None'    (Default) Return frequencies as counts.
         %
         %             o 'Rows'    Divide frequency counts by total
         %                         frequency in row (distribution within
         %                         row or within bin).
         %
         %             o 'Columns' Divide frequency counts by total
         %                         frequency in columns (distribution
         %                         within column or across bins).
         %
         %             o 'Total'   Divide frequency by total number of
         %                         observations in frequency table.
         %
         % Notes:
         %
         %   o This base-class implementation calls the method
         %     getBinLevelFrequencyTable that is implemented in the derived
         %     classes.
         
         p = inputParser;
         p.addParameter('BinLevel',true,@islogical);
         p.addParameter('Normalize','None',@(x)ismember(lower(x),...
            {'none','rows','columns','total'}));
         parse(p,varargin{:});
         BinLevel = p.Results.BinLevel;
         Normalize = p.Results.Normalize;
         
         if BinLevel
            [FreqTab, NanFreqTab] = obj.getBinLevelFrequencyTable;
         else
            FreqTab = obj.FrequencyTable;
            NanFreqTab = obj.nanFrequencyTable;
         end        
         
         if ~strcmpi(Normalize,'none')
            switch lower(Normalize)
               case 'rows'
                  FreqTab = bsxfun(@rdivide,FreqTab,sum(FreqTab,2));
               case 'columns'
                  FreqTab = bsxfun(@rdivide,FreqTab,sum(FreqTab,1));
               case 'total'
                  FreqTab = FreqTab./sum(FreqTab(:));
            end
         end
         
      end

      function Stats = getStatistics(obj,varargin)
         % GETSTATISTICS Compute statistics.
         %
         % Given a frequency table, there are many statistics that can be
         % computed. This method can compute statistics at bin level, or
         % observation level.
         % Optional Inputs:
         %
         % 'BinLevel'  Flag to indicate whether the statistics of interest
         %             are at bin level (true, default value), or at the
         %             observation level (if false).
         %
         % 'StatsList' Use this option to indicate which statistics need to
         %             be computed. Options are:
         %
         %             o 'All'    (Default) Return all supported
         %                         statistics.
         %
         %             o 'Odds'    For binary responses only, this is the
         %                         ratio of the first to the second column
         %                         of the frequency table.
         %
         %             o 'WOE'     For binary responses only, this is the
         %                         log of the distribution of observations
         %                         in column one, to the distribution of
         %                         observations in column two.
         %
         %             o 'InfoVal' For binary responses only, this is the
         %                         weighted sum of WOE values across bins,
         %                         where the weights are given by the
         %                         difference of the distribution in column
         %                         one minus the distribution in column two.
         %
         % Output:
         %
         %   o Stats: Structure with statistics stored in separate fields.
         %      For statistics that consist of a single, total amount such
         %      as information value, there is an extra field containing
         %      the contribution of each row to the total statistic.
         %
         
         p = inputParser;
         p.addParameter('BinLevel',true,@islogical);
         p.addParameter('StatsList',{'all'},@(x)all(ismember(lower(x),...
            {'all','infoval','odds','woe','entropy'})));
         parse(p,varargin{:});
         BinLevel = p.Results.BinLevel;
         if ischar(p.Results.StatsList)
            StatsList = {p.Results.StatsList};
         else
            StatsList = p.Results.StatsList;
         end
         StatsList = lower(StatsList(:));
         
         [FreqTab, nanFreqTab] = obj.getFrequencyTable('BinLevel',BinLevel);
         Stats = GetStatisticsByFreq_mod(FreqTab,nanFreqTab,...
            'StatsList',StatsList);
         
      end
      
      function obj = setLabels(obj,labels)
         % SETLABELS Set bin labels.
         %
         % This base-class implementation handles labels setting for
         % derived classes.
         
         if ~iscell(labels) || ~all(cellfun(@(x)ischar(x),labels))...
               || ~isvector(labels)
            error(message('finance:internal:finance:binning:Container:InvalidLabelsType'));
         end
         
         if length(labels)~=obj.NumBins
            error(message('finance:internal:finance:binning:Container:InvalidLabelsLength'));
         end
         
         if numel(unique(labels(:)))~=numel(labels(:))
            error(message('finance:internal:finance:binning:Container:DuplicateLabels'));
         end
         
         obj.UserLabels = labels(:);
         
      end
      
   end
   
   methods (Abstract = true)
      
      labels = getLabels(obj) % Get method for labels, implemented in derived classes

   end
   
   methods (Access = protected, Abstract = true)
      
      ft = getBinLevelFrequencyTable(obj) % Get bin frequencies, implemented in derived classes

   end
   
end
