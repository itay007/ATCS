function s = featuresparse(featxt,varargin)
% FEATURESPARSE parses features from GenBank, GenPept, or, EMBL data.
%
%   FTRS = FEATURESPARSE(STR) parses the features from the input string STR
%   which contains a GenBank, GenPept, or EMBL section of features. STR may
%   also be a char array with the text respective to the features or a
%   structure such as the returned by GENBANKREAD, GENPEPTREAD, or,
%   EMBLREAD. FTRS is the output structure, where every field name 
%   corresponds to the name of every feature. Fields contain sub-structures
%   with feature qualifiers as fields. The only mandatory qualifier is the
%   'Location'. When possible, 'Location' is also translated to numeric
%   indices and stored into the field 'Indices'. When using the numeric
%   indices to extract sequence information observe that you may still need
%   to complement the sequences.
%
%   Fields in the returning structure will have the same name as in the
%   sequence database all except the following special cases:
%
%   -10_signal  =>  minus_10_signal
%   -35_signal  =>  minus_35_signal
%   3'UTR       =>  three_prime_UTR
%   3'clip      =>  three_prime_clip
%   5'UTR       =>  five_prime_UTR
%   5'clip      =>  five_prime_clip
%   D-loop      =>  D_loop
%
%   FTR = FEATURESPARSE(...,'Feature',FEATURE_NAME) returns only the
%   sub-structure respective the FEATURE_NAME. If there are several 
%   features with the same name FTR is an array of structures. 
%
%   FTR = FEATURESPARSE(...,'Sequence',true) extracts, when possible, the
%   sequences respective to each feature, joining and complementing pieces
%   of the source sequence and storing them into the field 'Sequence'. When
%   extracting incomplete CDS, FEATURESPARSE uses the codon_start qualifier
%   to adjust the frame of the sequence. Default is false. 
%
%   Examples:
% 
%       % Obtain all the features stored in the file 'nm175642.txt':
%       gbk_struct = genbankread('nm175642.txt');  
%       features = featuresparse(gbk_struct)
%
%       % Obtain the CDS of the Caenorhabditis elegans cosmid record in
%       % genbank:
%       worm = getgenbank('Z92777');
%       CDS = featuresparse(worm,'feature','cds')
%
%   See also EMBLREAD, GENPEPTREAD, GENBANKREAD, GETGENBANK, GETGENPEPT,
%   SEQSTATSDEMO.

% Copyright 2003-2012 The MathWorks, Inc.


%   Hidden option:
%   FEATURESPARSE(...,'Text',true) adds a field named 'Text' to the
%   sub-structures with the text of every feature.

% default options
saveText = false;
onlyOneFeature = false;
trimNans = false;
getSequences = false;
sourceSequence = [];

% process input arguments
if  nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:featuresparse:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'feature','text','trimnans','sequence'};
    for j=1:2:nargin-1
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname))); 
        if isempty(k)
            error(message('bioinfo:featuresparse:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:featuresparse:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % 'feature'
                    onlyOneFeature = true;
                    requiredFeature = pval;
                    if ~ischar(pval)
                        error(message('bioinfo:featuresparse:InvalidString'));
                    end
                case 2 % 'text'
                    saveText = bioinfoprivate.opttf(pval);
                    if isempty(saveText)
                        error(message('bioinfo:featuresparse:verboseInputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 3 % 'trimnans'
                    trimNans = bioinfoprivate.opttf(pval);
                    if isempty(trimNans)
                        error(message('bioinfo:featuresparse:verboseInputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 4 % 'sequence'
                    getSequences = bioinfoprivate.opttf(pval);
                    if isempty(getSequences)
                        error(message('bioinfo:featuresparse:verboseInputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
            end
        end
    end
end

% Check if the input is a structure which has the field 'features' (as in
% genbankread and genpeptread)
if isfield(featxt,'Features')
    if numel(featxt)>1
        error(message('bioinfo:featuresparse:NoArrayStructuresGenbank'))
    end
    if getSequences
        sourceSequence = featxt.Sequence;
    end
    featxt = featxt.Features;
end
% Check if the input is a structure which has the field 'Feature' (as in
% emblread)
if isfield(featxt,'Feature') 
    if numel(featxt)>1
        error(message('bioinfo:featuresparse:NoArrayStructuresEmblread'))
    end
    if getSequences
        sourceSequence = featxt.Sequence;
    end
    featxt = featxt.Feature;
end

if getSequences && isempty(sourceSequence)
    error(message('bioinfo:featuresparse:noSourceSequence'))
end

% Check if the input is a single string with linefeed characters
if ischar(featxt) && size(featxt,1)==1
    featxt = char(strread(featxt,'%s','delimiter','\n','whitespace',''));
end

% Now feattxt can only be a char array
if ~ischar(featxt)
    error(message('bioinfo:featuresparse:InvalidInputArgument'))
end

%initialize output structure
s = struct;

%line counter
ln = 1;

% check if the very first line is the header of the feature section
if ~isempty(regexp(featxt(ln,:),'Location/Qualifiers','once'))
    ln = ln+1;
end

moreFeatures = true;
while moreFeatures % main loop, searches for features
    
    % find name and location of next feature
    feaNameId = regexp(featxt(ln,:),'\S+','once');
    feaName   = regexp(featxt(ln,:),'\S+','match','once');
    
    % special case for 5'UTR and 3'UTR (field names cannot start with
    % number)
    if ~isempty(regexp(feaName,'^\d''UTR$','once')) 
        switch feaName(1)
            case '3'
                feaName = 'three_prime_UTR';
            case '5'
                feaName = 'five_prime_UTR';
        end
    end
    % special case for 5'clip and 3'clip (field names cannot start with
    % number)
    if ~isempty(regexp(feaName,'^\d''clip$','once')) 
        switch feaName(1)
            case '3'
                feaName = 'three_prime_clip';
            case '5'
                feaName = 'five_prime_clip';
        end
    end
    % special case for -10_signal and -35_signal (field names cannot start with
    % number)
    if ~isempty(regexp(feaName,'^-\d\d_signal$','once')) 
        feaName = ['minus_' feaName(2:10)];
    end
    
    % special case for placeholder (field names cannot start with -)
    if ~isempty(regexp(feaName,'^-$','once')) 
        feaName = 'located_sequence_feature';
    end
    
    % special case for placeholder (field names cannot have a -)
    if ~isempty(regexp(feaName,'-','once')) 
        feaName = regexprep(feaName,'-','_');
    end
        
    % find if we already have a field with that feature name
    if isfield(s,feaName)
        feaCounter = numel(s.(feaName))+1;
    else
        feaCounter = 1;
    end
    
    % start of this feature (for saving later the whole text)
    se_ln = ln;
    
    % The rest of the line should be the location
    feaLocation = featxt(ln,(feaNameId+numel(feaName)+1):end);
    feaLocation = regexprep(feaLocation,'\s','');
    ln = ln + 1;
    
    % using parentheses to find if there are more lines belonging to the
    % location 
    while sum(feaLocation=='(') > sum(feaLocation==')')
        feaLocation = [feaLocation regexprep(featxt(ln,:),'\s','')]; %#ok<AGROW>
        ln = ln + 1;    
    end
    
    if ~onlyOneFeature||(onlyOneFeature && strcmpi(requiredFeature,feaName))
        % store location into output structure
        s.(feaName)(feaCounter).Location = feaLocation;
        [Indices,UFB] = featurelocation(feaLocation);
        if trimNans
            s.(feaName)(feaCounter).Indices = Indices(~isnan(Indices));
        else
            s.(feaName)(feaCounter).Indices = Indices;
        end
        if UFB
            s.(feaName)(feaCounter).UnknownFeatureBoundaries = true;
        end
    end
    
    while true % inner loop: searches for feature qualifiers
        % qualifier names must start with /, check if we have a qualifier,
        % we may also be at the end of the text of at the beginning of other
        % feature
        if ln<=size(featxt,1)
            qualifierNameId = regexp(featxt(ln,:),'(?<=/)[^=]+','once');
        else
            qualifierNameId = [];
        end
        if isempty(qualifierNameId)
            break;
        end
        
        % read in the qualifier name
        qualifierName = regexp(featxt(ln,qualifierNameId:end),'[^=\s]+','match','once');
        
        % get some tokens of the rest of the line to figure out the type of qualifier value
        t = regexp(featxt(ln,qualifierNameId+numel(qualifierName):end),'^\s*(=?)\s*("?)(.*)$','tokens','once');
        % t{1} gets the = (if any)
        % t{2} gets the first " (if any)
        % t{3} gets the rest of the line
        
        if isempty(t{1}) % CASE 1: catch the case with no qualifier value: /qualifier_name
            %if ~isempty(t{2}) || ~isempty(t{3})
                %do not error or warn just continue silently
            %end
            qualifierValue = true;
            ln = ln + 1;
        elseif isempty(t{2}) % CASE 2:  qualifier value is not free text (single or multi line)
            qualifierValue = strtrim(t{3});
            ln = ln+1;
            % using parentheses to find if there are more lines belonging to this qualifier value
            while sum(qualifierValue=='(') > sum(qualifierValue==')')
                qualifierValue = [qualifierValue regexprep(featxt(ln,:),'\s','')]; %#ok<AGROW>
                ln = ln + 1;
            end           
        else % CASE 3: qualifier value is free text surrounded by "" (single or multi line)
            qualifierValueRow = {};
            thisRow = strtrim(t{3});
            rowNumber = 1;
            ln = ln+1;
            % line with the end of the qualifier must have ",""","""""",...
            % preventing the case where a escaped " is at the end of the line
            while  ~rem(numel(regexp(thisRow,'"*(?=\s*$)','match','once')),2)
                qualifierValueRow{rowNumber} = thisRow;   %#ok<AGROW>
                thisRow = strtrim(featxt(ln,:));
                rowNumber = rowNumber+1;
                ln = ln+1;
            end
            qualifierValueRow{rowNumber} = regexprep(thisRow,'"\s*$',''); %#ok<AGROW>
            % should we put spaces within every row ? only if we find a
            % space in any of the rows otherwise we concatenate without
            % extra  spaces
            if sum(cellfun('length',regexp(qualifierValueRow,'\s')))
                qualifierValue = strtrim(sprintf('%s ',qualifierValueRow{:}));
            else
                qualifierValue = sprintf('%s',qualifierValueRow{:});
            end            
            % replace all escaped "
            qualifierValue = regexprep(strtrim(qualifierValue),'""','"');
        end
        
        if ~onlyOneFeature||(onlyOneFeature && strcmpi(requiredFeature,feaName))
            
            % now store the qualifier info in to the structure:
            % does a field with the qualifier name already exist and is not empty ?
            % (it is possible that the field was created because features
            % are stored as array of structures)
            if isfield(s.(feaName)(feaCounter),qualifierName) && ~isempty(s.(feaName)(feaCounter).(qualifierName))
                % if it exist and it is a cell then it is already a repeated qualifier
                if iscell(s.(feaName)(feaCounter).(qualifierName))
                    qualifierCounter = numel(s.(feaName)(feaCounter).(qualifierName))+1;
                else % if it is not a cell then convert it to a cell to store several values
                    s.(feaName)(feaCounter).(qualifierName) = {s.(feaName)(feaCounter).(qualifierName)};
                    qualifierCounter = 2;
                end
            else
                qualifierCounter = 1;
            end
            
            switch qualifierCounter
                case 1
                    s.(feaName)(feaCounter).(qualifierName)=qualifierValue;
                otherwise
                    s.(feaName)(feaCounter).(qualifierName){qualifierCounter}=qualifierValue;
            end
            
        end
        
        
    end % while more qualifiers

    % extract the sequence when necessary
    if ~onlyOneFeature||(onlyOneFeature && strcmpi(requiredFeature,feaName))
        % Indices and UFB are still valid here
        if getSequences && ~any(isnan(Indices)) && ~isempty(Indices)
            numSegments = numel(Indices)/2;
            Indices = reshape(Indices,2,numSegments);
            seq = blanks(sum(abs(diff(Indices,[],1))+1));
            j = 1;
            for i = 1:numSegments
                if Indices(1,i)<=Indices(2,i)
                    k = j+Indices(2,i)-Indices(1,i);
                    seq(j:k) = sourceSequence(Indices(1,i):Indices(2,i));
                else
                    k = j+Indices(1,i)-Indices(2,i);
                    seq(j:k) = seqrcomplement(sourceSequence(Indices(2,i):Indices(1,i)));
                end
                j = k+1;
            end
            %special case for incomplete coding regions that start out of frame
            if UFB && strcmp(feaName,'CDS') && isfield(s.CDS(feaCounter),'codon_start')
                seq = seq(str2double(s.CDS(feaCounter).codon_start):end);
            end
            s.(feaName)(feaCounter).Sequence = seq;
        end
    end
    
    
    % save the text of the features into structure
    if saveText
        s.(feaName)(feaCounter).Text = featxt(se_ln:ln-1,:);
    end
    
    % end of features text?
    if ln>size(featxt,1)
        moreFeatures = false;
    end
end

if onlyOneFeature
   fieldNames = fieldnames(s);
   h = find(strcmpi(fieldNames,requiredFeature));
   if isempty(h)
       s = struct('');
       warning(message('bioinfo:featuresparse:featureNotFound', requiredFeature))
   else
       s = s.(fieldNames{h(1)});
   end
end



