function I = seqmatch(strs,lib,varargin)
%SEQMATCH Find matches for every string in a library.
%   IND = SEQMATCH(STRS,LIB) looks through the elements of LIB to find
%   strings that begin with every string in STRS. IND contains the indices
%   the first occurrence for every string in the query. If no match is found
%   for a given query the respective index is 0. STRS and LIB must be cell
%   arrays of strings. 
%
%   IND = SEQMATCH(STRS,LIB,'exact',true) looks for exact matches only.
%
%   Example:
%      
%     lib = {'VIPS_HUMAN','SCCR_RABIT','CALR_PIG','VIPR_RAT','PACR_MOUSE'};
%     query = {'CALR','VIP'};
%     h = seqmatch(query,lib);
%     lib(h)
%
%   See also REGEXPI, SEQSHOWWORDS, SEQWORDCOUNT.

%   Copyright 2003-2012 The MathWorks, Inc.


doExactMatch = false;

if  nargin > 2
    if rem(nargin,2) == 1
        error(message('bioinfo:seqmatch:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'exact',''};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname))); 
        if isempty(k)
            error(message('bioinfo:seqmatch:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:seqmatch:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1
                   doExactMatch = bioinfoprivate.opttf(pval);
                    if isempty(doExactMatch)
                        error(message('bioinfo:seqmatch:InputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
            end %switch
        end %if
    end %for
end %if
         
if ~iscell(strs)
    error(message('bioinfo:seqmatch:IncorrectSTRSType'))
end

if ~iscell(lib)
    error(message('bioinfo:seqmatch:IncorrectLIBType'))
end

strs = strs(:);
lib = lib(:);

numStrs = numel(strs);

if ~iscellstr(strs)
    error(message('bioinfo:seqmatch:IncorrectSTRSElementType'))
end
if ~iscellstr(lib)
    error(message('bioinfo:seqmatch:IncorrectLIBElementType'))
end

I = zeros(numStrs,1);
for i = 1:numStrs
   if doExactMatch
       h = strmatch(strs{i},lib,'exact'); 
   else
       h = strmatch(strs{i},lib); 
   end
   if ~isempty(h)
       I(i) = h(1);
   end
end

if any(I==0)
     warning(message('bioinfo:seqmatch:StringNotFound'))
end

 
