classdef Algorithm_mod
% ALGORITHM Base class for automatic binning algorithms.
%
%   Automatic binning algorithms can be implemented as derived classes from
%   Algorithm. The Algorithm class handles the algorithm options in a
%   generic way to simplify the implementation of particular algorithms.
%
%   Derived classes use addOption to add options for the specific algorithm
%   they implement. The Algorithm base class can then handle options
%   setting (setOptions), getting (getOptions) and resetting to defaults
%   (resetDefaultOptions).
%
%   The runAlgorithm method is abstract, and that is where derived
%   classes actually implement a specific automatic binning algorithm.
%
%   This is an internal function that is not meant to be called directly
%   by the user.

% Copyright 2014 The MathWorks, Inc.
   
   properties (SetAccess = protected)
      AlgorithmID = ''; % Algorithm identifier
      IsSupervised = []; % The algorithm supervised (true) or unsupervised (false)
   end

   properties (Access = protected)
      Options = struct; % Current algorithm options
      DefaultOptions = struct; % Default algorithm options
      ValidationFunctions = struct; % Options' validation functions
   end

   methods
      
      function obj = Algorithm_mod(algoID)
         % ALGORITHM Create an Algorithm object.
         obj.AlgorithmID = algoID;
      end
      
      function obj = setOptions(obj,varargin)
         % SETOPTIONS Set algorithm options.
         %
         % Options come in name-value pairs. The names must match the
         % option names the Algorithm knows about (Options property). The
         % default values come from the default options stored in Algorithm
         % (DefaultOptions property). The validation of option values is
         % performed using the validation functions stored in Algorithm
         % (ValidationFunctions property).
         %
         % This base-class implementation handles option setting for
         % derived classes.
         %
         % The option names, default values and validation functions are
         % set by derived classes with the addOption method.
         
         optnames = fields(obj.Options);
         
         if isempty(optnames)
            warning(message('finance:internal:finance:binning:Algorithm:NoOptions'));
            return;
         end
         
         p = inputParser;
         
         for ii=1:length(optnames)
            name = optnames{ii};
            % By setting current option as default value in parser, all
            % options previously set are honored
            val = obj.Options.(name);
            valfn = obj.ValidationFunctions.(name);
            p.addParameter(name,val,valfn);
         end
         
         try
             p.parse(varargin{:});
         catch MException
             throw(MException)
         end
         
         for ii=1:length(optnames)
            name = optnames{ii};
            obj.Options.(name) = p.Results.(name);
         end
         
      end
      
      function opt = getOptions(obj,optname)
         % GETOPTIONS Get algorithm options.
         %
         % Returns a structure with the current Algorithm options.
         %
         % This base-class implementation handles option getting for
         % derived classes.
         
         if nargin==1
            opt = obj.Options;
            return;
         end
         
         optnames = fields(obj.Options);
         if isempty(optnames)
            warning(message('finance:internal:finance:binning:Algorithm:NoOptions'));
            return;
         end
         
         [ok,ind] = ismember(lower(optname),lower(optnames));
         if ok
            name = optnames{ind};
            opt = obj.Options.(name);
         end

      end

      function obj = resetDefaultOptions(obj)
         % RESETDEFAULTOPTIONS Reset options to default values.
         %
         % This base-class implementation handles resetting to default
         % options for derived classes.
         
         optnames = fields(obj.Options);
         if isempty(optnames)
            warning(message('finance:internal:finance:binning:Algorithm:NoOptions'));
            return;
         end
         
         for ii=1:length(optnames)
            name = optnames{ii};
            obj = obj.setOptions(name,obj.Options.(name));
         end

      end
      
      function displayOptions(obj)
         % DISPLAYOPTIONS Display options in command line.
                  
         fprintf('\nAlgorithm options ');
         fprintf(' ( Option Name: Value (Default) )\n\n');
         optnames = fields(obj.Options);
         for ii=1:length(optnames)
            name = optnames{ii};
            val = obj.Options.(name);
            defval = obj.DefaultOptions.(name);
            if ischar(val)
               fprintf('   %s: %s (%s)\n',name,val,defval);
            else
               fprintf('   %s: %g (%g)\n',name,val,defval);
            end
         end
         fprintf('\n');
         
      end
      
      function dbc = applySorting(~,dbc,Sorting)
          % Apply sorting 'Sorting' to the categories in the binning
          % container 'dbc'. 'dbc' must be of type Categorical/Nominal.
          % This method is called from Algorithm subclasses. The values of
          % the input parameter 'Sorting' are:
          %   'Odds'   : The category groupings, and thus the frequency 
          %              table are sorted by order increasing values of 
          %              odds.
          %   'Goods'  : The category groupings, and thus the frequency 
          %              table are sorted by order increasing values of 
          %              Goods. The indexing is mapped to the column of
          %              Bads in the frequency table.
          %   'Bads'   : The category groupings, and thus the frequency 
          %              table are sorted by order increasing values of 
          %              Bads. The indexing is mapped to the column of
          %              Goods in the frequency table.
          %   'Totals' : The category groupings, and thus the frequency 
          %              table are sorted by order increasing values of 
          %              Goods. The indexing is mapped to the column of
          %              Bads in the frequency table.
          %   'None'   : No sorting is applied.
          
          switch lower(Sorting)
              case 'odds'
                  dbc = dbc.sortPredictorCategories('Criterion',...
                      'odds','Mode','ascending');
              case 'goods'
                  dbc  = dbc.sortPredictorCategories('Criterion',...
                      'colcount','ColNumber',1,'Mode','ascending');
              case 'bads'
                  dbc  = dbc.sortPredictorCategories('Criterion',...
                      'colcount','ColNumber',2,'Mode','ascending');
              case 'totals'
                  dbc  = dbc.sortPredictorCategories('Criterion','RowCount');
          end
      end
      
   end

   methods (Access = protected)

      function obj = addOption(obj,name,defval,validfn)
         % ADDOPTION Add an algorithm option.
         %
         % Derived classes set option names, default values and validation
         % functions by calling this method.
         %
         % This method should be called at construction of the derived
         % class, to indicate the algorithm's option information to the
         % base class, and from that point the base class is responsible
         % for setting, getting and resetting options.
         %
         % Note that the access is protected, only derived classes can add
         % options, clients of particular instances of an algorithm cannot
         % add options.
         
         obj.Options.(name) = defval;
         obj.DefaultOptions.(name) = defval;
         obj.ValidationFunctions.(name) = validfn;
         
      end
      
      function dbc = setBinsByCutPointIndices(~,dbc,CutPointIndices, nanBin)
         % SETBINSBYCUTPOINTINDICES Use cut-point _indices_ to modify the
         %   bins of a data binning container object.
         %
         % Algorithms often keep track of bin boundaries with pointers to
         % the rows of the data's frequency table. The term 'cut-point
         % indices' refers to these pointers. However, the data binning
         % container interface requires setting bin boundaries as cut
         % points (not there indices) for numeric data, and as category
         % groupings (not a mapping to row indices) for categorical
         % information.
         %
         % This method is a utility to transform cut-point indices to cut
         % points for numeric data, or to category groupings for
         % categorical data. It allows derived algorithms to simply call
         % this method to update the bins of a container, instead of having
         % to implement the conversion to cut points or category groupings
         % themselves.
         
         if isa(dbc,'Numeric_mod')
            if isempty(CutPointIndices)
               dbc = dbc.setCutPoints([], nanBin);
            else
               CutPoints = dbc.SortedPredictorValues(CutPointIndices);
               dbc = dbc.setCutPoints(CutPoints, nanBin);
            end
         else
            CatGrouping = dbc.getCatGrouping;
            NumRows = size(CatGrouping,1);
            BinNumber = ones(NumRows,1);
            if ~isempty(CutPointIndices)
               NumCuts = length(CutPointIndices);
               for ii = 2:NumCuts
                  i1 = CutPointIndices(ii-1);
                  i2 = CutPointIndices(ii)-1;
                  BinNumber(i1:i2) = ii;
               end
               BinNumber(CutPointIndices(end):end) = NumCuts+1;
            end
            CatGrouping.BinNumber = BinNumber;
            dbc = dbc.setCatGrouping(CatGrouping);
         end

      end
      
   end
   
   methods (Abstract = true)
      
      dbc = runAlgorithm(obj,dbc) % Implementation of the actual algorithm
      
   end
   
end
