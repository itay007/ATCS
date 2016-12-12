function varargout = browserFunc(whichcall, varargin)
%BROWSERFUNC  Support function for the NGSBrowser component

%   Copyright 2010-2012 The MathWorks, Inc.

persistent genomeBrowserDataMap; 
if ~isa(genomeBrowserDataMap,'containers.Map')
    genomeBrowserDataMap = containers.Map;
end

varargout = cell(1, 1);
switch whichcall
    case 'validateAndLoadBioMapObjectFromWS'
        % called from: GenomeController.addDataTrackFromMatlab
        [genomeBrowserDataMap, varargout{1},varargout{2}] = validateAndLoadBioMapObjectFromWS(genomeBrowserDataMap, varargin{1});
    case 'getBioDataVariables'
        % called from: LoadBioObjectDataDialog.getVariablesFromMatlab
        [varargout{1}, varargout{2}] = getBioDataVariables();
    case 'clearAll'
        % called from: GenomeController.clearAllBioDataMap
        genomeBrowserDataMap = [];
    case 'getSequenceDictionary' 
        % called from: GenomeMenuBar.addShortReadDataFromFile
        [varargout{1}] = getSequenceDictionary(varargin{1});
    case 'getScannedDictionary'
        % called from: SelectReferenceDialog.getScannedDictionaryFromMatlab
        [varargout{1}, varargout{2}, varargout{3}] = getScannedDictionary(varargin{1});        
    case 'getDataRange'
        % called from: AbstractTrack.loadDataRange
        %              AlignmentTrack.loadDataRange
        [genomeBrowserDataMap, varargout{1},varargout{2}, varargout{3}, varargout{4}, varargout{5}] = getDataRange(genomeBrowserDataMap, varargin(:));
    case 'clearBioDataObject'
        % called from: AbstractTrack.unload
        genomeBrowserDataMap = clearBioDataObject(genomeBrowserDataMap, varargin{1});        
    case 'getCompactAlignment'
        % called from: AlignmentTrack.getAlignmentFromMatlab
        [varargout{1}, varargout{2}, varargout{3}] = getPreProcessedAlignments(genomeBrowserDataMap, varargin{1}, varargin{2}, varargin{3});
    case 'getFastaSequenceInfo'
        [varargout{1}, varargout{2}] = getFastaSequenceInfo(genomeBrowserDataMap, varargin{1});
    case 'getFastaSequenceAt',
        % called from: SequenceTrack.getSequenceBlockFromMatlab
        [varargout{1},varargout{2},varargout{3}, varargout{4}] = getFastaSequenceAt(genomeBrowserDataMap, varargin{1}, varargin{2}, varargin{3});
    case 'getCoverage'
        %called from: AlignmentTrack.getCoverageFromMatlab
        [varargout{1}, varargout{2}, varargout{3}, varargout{4}, varargout{5}] = reportCoverage(genomeBrowserDataMap, varargin{1}, varargin{2}, varargin{3}, varargin{4});
    case 'getFeatures'
        %called from: FeaturesTrack.loadFeatureFromMatlab
        [varargout{1}, varargout{2}, varargout{3}, varargout{4}] = getFeatures(genomeBrowserDataMap, varargin{1}, varargin{2}, varargin{3});
    otherwise
        error(message('bioinfo:browserFunc:UnknownOption', upper( whichcall )));
end
end

function [genomeBrowserDataMap,bioMapKey,errMsg] = validateAndLoadBioMapObjectFromWS(genomeBrowserDataMap,bioMapName)
    
    n = numel(bioMapName);
    bioMapKey = repmat({''},1,n);
    errMsg = repmat({''},1,n);
    
    for i = 1:n
        bioMapVar = evalin('base', bioMapName{i});
    
        % Validate BioMap object:
        if ~(isobject(bioMapVar) && isa(bioMapVar, 'BioMap'))
            errMsg{i} = getString(message('bioinfo:ngsbrowser:NotAValidBioMapObject'));
            continue;
        elseif numel(bioMapVar.SequenceDictionary)>1
            errMsg{i} = getString(message('bioinfo:ngsbrowser:BioMapDataMultipleReferences',bioMapName{i}));
            continue;
        elseif bioMapVar.NSeqs==0
            errMsg{i} = getString(message('bioinfo:ngsbrowser:BioMapDataEmpty',bioMapName{i}));
            continue;
        else
            st = getStart(bioMapVar);
            if ~issorted(st(st>0))
                errMsg{i} = getString(message('bioinfo:ngsbrowser:BioMapDataNotSorted',bioMapName{i}));
                continue;
            end
        end    
        bioMapKey{i} = getGenomeDataMapKey(genomeBrowserDataMap, bioMapName{i});
        genomeBrowserDataMap(bioMapKey{i}) = bioMapVar;
        
    end
    
end

function [genomeBrowserDataMap, rangeStart,rangeEnd,errMsg, mapKey, isMapEmpty] = getDataRange(genomeBrowserDataMap, varargin)
errMsg = '';
isMapEmpty  = false;
isReference = false;
mapKey = varargin{:}{1};
[~, ~, extStr] = fileparts(mapKey);
try
    if exist(mapKey, 'file') && ~isempty(extStr);
        %Check if it is a reference file string
        dataObj = [];
        if ~isempty(regexp(lower(mapKey),'\.fasta$','once')) || ...
           ~isempty(regexp(lower(mapKey),'\.fa$','once')) || ...
           ~isempty(regexp(lower(mapKey),'\.fas$','once')) || ...
           ~isempty(regexp(lower(mapKey),'\.fna$','once')) 
            dataObj =  bioinfoprivate.MemoryMappedFastaFile(mapKey);
            rangeStart =  dataObj.Range(1);
            rangeEnd = dataObj.Range(2);
            isReference = true;
        elseif ~isempty(regexp(lower(mapKey),'\.gff$','once'))
            dataObj = GFFAnnotation(mapKey);
            dataRange = dataObj.getRange;
            rangeStart =  dataRange(1);
            rangeEnd = dataRange(2);
        elseif ~isempty(regexp(lower(mapKey),'\.gtf$','once'))
            dataObj = GTFAnnotation(mapKey);
            dataRange = dataObj.getRange;
            rangeStart =  dataRange(1);
            rangeEnd = dataRange(2);
        elseif ~isempty(regexp(lower(mapKey),'\.sam$','once'))
            if numel(varargin{:})==1
               dataObj = BioMap(mapKey);
            else 
                dataObj = BioMap(mapKey, 'SelectReference', varargin{:}{2});
            end
            if dataObj.NSeqs == 0
                % Empty BioMap
                isMapEmpty = true;
                rangeStart = 0;
                rangeEnd = 0;
                return;
            end               
            [rangeStart,rangeEnd] = dataObj.getRange;
        elseif ~isempty(regexp(lower(mapKey),'\.bam$','once'))
            if numel(varargin{:})==1
                dataObj = BioMap(mapKey);
            else
                dataObj = BioMap(mapKey, 'SelectReference', varargin{:}{2});
            end
            if dataObj.NSeqs == 0
                % Empty BioMap
                isMapEmpty = true;
                rangeStart = 0;
                rangeEnd = 0;
                return;
            end              

            [rangeStart,rangeEnd] = dataObj.getRange;
        end
        
        if ~isempty(dataObj)
            [~, fileName]= fileparts(mapKey);
            if ~isReference
                fileName = getGenomeDataMapKey(genomeBrowserDataMap, fileName);
            else
                fileName = ['Reference-', fileName];
            end
            genomeBrowserDataMap(fileName) = dataObj;
            mapKey = fileName;
        else
            error(message('bioinfo:browserFunc:UnknowFileType', mapKey));
        end
    else
        obj= genomeBrowserDataMap.values({mapKey});
        obj = obj{:};
        [rangeStart,rangeEnd] = obj.getRange;
    end
catch ME
    errMsg = ME.message;
    rangeStart = 0;
    rangeEnd=0;
end
end

function [alignments, order, err] = getPreProcessedAlignments(genomeBrowserDataMap, bioMapKey, startPos, endPos)
   err = '';
   alignments = [];
   order = [];
   try
        obj= genomeBrowserDataMap.values({bioMapKey});
        obj = obj{:};
        startPos = startPos + 1;
        endPos = endPos + 1;
        % Check not exceed the limits
        [minStart, maxEnd] = obj.getRange;
        
        if startPos < minStart
            startPos = minStart;
        end
        if startPos == 0
            startPos = 1;
        end
        
        if endPos > maxEnd
            endPos = maxEnd;
        end
        
        [order, idx] = getRowInCompactAlignment(obj, double(startPos), double(endPos), false);
        
        if ~isempty(order)
            % Get alignments Java array
            alignments = alignmentRecords(getSubset(obj, idx));
            % Find number of unique orders and add it to the top of the order list
            order = order(:);
            order = [numel(unique(order)); order];
        else
            err = 'No reads found';
        end 
   catch ME
       err = ME.message;
   end
end

function alignments = alignmentRecords(bmObj)
% Return a javaArray of alignments. Each alignment is a Java object
% Alignment containing read information for each alignment in the
% alignIndices.
import com.mathworks.toolbox.bioinfo.genome.data.*;
% Get properties about the records
data = get(bmObj,{'Signature','Quality','Header','Sequence'});
Headers = data{3};
Signatures = data{1};
GappedSequences = bioinfoprivate.cigar2gappedsequence(data{4}, Signatures, 'GapsInRef', true);
Starts = int32(getStart(bmObj));
Ends = int32(getStop(bmObj));
MappingQualities = int8(getMappingQuality(bmObj));
if ~isempty(data{2})
    GappedQualities = bioinfoprivate.cigar2gappedsequence(data{2}, Signatures, 'quality', true);
    % Base quality as ASCII, convert tem to binary phred scores
    GappedQualities = cellfun(@double,GappedQualities,'UniformOutput', false);
else
    GappedQualities = [];
end
% Flags
readNegativeStrandFlags = filterByFlag(bmObj, 'strandQuery', 1); % 1 is reverse
duplicateFlags = filterByFlag(bmObj, 'duplicate', true);
readUnmappedFlags = filterByFlag(bmObj, 'unmappedQuery', true);
readPairedFlags = filterByFlag(bmObj, 'pairedInSeq', true);
mateNegativeStrandFlags = filterByFlag(bmObj,'strandMate', 1);
mateUnmappedFlags = filterByFlag(bmObj, 'unmappedMate', true);
properPairFlags = filterByFlag(bmObj, 'pairedInMap', true);

% Pair-mate info
mateAlignmentStarts = zeros(1,bmObj.NSeqs); %% Need to get
mateReferenceNames = Headers; %% Need to get

% Construct javaArray
alignments = javaArray('com.mathworks.toolbox.bioinfo.genome.data.ReadAlignment', bmObj.NSeqs);
if isempty(GappedQualities)
    for i = 1:bmObj.NSeqs
        alignments(i) = ReadAlignment(['test_' num2str(i)], Starts(i), Ends(i), ...
            10, GappedSequences{i}, Signatures{i},...
            false, false,...
            false, false);
    end
    
else
    for i = 1:bmObj.NSeqs
        if readPairedFlags(i)
            alignments(i) = ReadAlignment(Headers{i}, Starts(i), Ends(i), ...
                MappingQualities(i), GappedSequences{i}, Signatures{i},...
                GappedQualities{i,:},...
                readNegativeStrandFlags(i), duplicateFlags(i), ...
                readUnmappedFlags(i), readPairedFlags(i), ...
                properPairFlags(i), mateReferenceNames{i}, mateAlignmentStarts(i), ...
                mateNegativeStrandFlags(i), mateUnmappedFlags(i));
        else %% It does not have a mate
            alignments(i) = ReadAlignment(Headers{i}, Starts(i), Ends(i), ...
                MappingQualities(i), GappedSequences{i}, Signatures{i},...
                GappedQualities{i,:},...
                readNegativeStrandFlags(i), duplicateFlags(i),...
                readUnmappedFlags(i), readPairedFlags(i));
        end
    end
end

end


function [coverage, startPos, endPos, maxCount, errMsg] = reportCoverage(genomeBrowserDataMap, varargin)
errMsg = '';
try
    bioMapKey = varargin{1};
    startPos = varargin{2};
    endPos = varargin{3};
    binWidth = varargin{4};
    
    obj= genomeBrowserDataMap.values({bioMapKey});
    obj = obj{:};
    
    [minStart, maxEnd] = obj.getRange;
    
    if startPos < minStart
        startPos = minStart;
    end
    
    if endPos > maxEnd
        endPos = maxEnd;
    end
    if binWidth <= 1
        coverage = getBaseCoverage(obj,  double(startPos), double(endPos));
    else
        coverage = getBaseCoverage(obj, double(startPos), double(endPos), 'BinWidth', binWidth);
    end
    maxCount = max(coverage);
catch ME
    errMsg = ME.message;
    coverage = [];
    maxCount = 0;
    startPos = 0;
    endPos = 0;
end
end

function [features, maxRow, maxGroup, errMsg] = getFeatures(genomeBrowserDataMap, varargin)
errMsg = '';
maxRow = 0;
maxGroup = 0;
try
    featureKey = varargin{1};
    startPos = varargin{2};
    endPos = varargin{3};
    
    obj= genomeBrowserDataMap.values({featureKey});
    obj = obj{:};
    
    dataRange = obj.getRange;
    
    if startPos < dataRange(1)
        startPos = dataRange(1);
    end
    
    if endPos > dataRange(2)
        endPos = dataRange(2);
    end
    
    features = obj.getBrowserData(double(startPos), double(endPos));
    if(~isempty(features))
        maxRow = max([features.RowInView]);
        maxGroup = max([features.GroupIndex]);
    end
catch ME
    errMsg = ME.message;
    features = [];
    maxRow = 0;
    maxGroup = 0;
end
end


function [refLength, refHeader] = getFastaSequenceInfo(genomeBrowserDataMap, refSeqKey)  
try
    obj= genomeBrowserDataMap.values({refSeqKey});
    obj = obj{:};
    if isstruct(obj)    
        refLength = numel(obj.Sequence);
        refHeader = obj.Header;
    else
        refLength = obj.Length;
        refHeader = obj.Header;
    end
catch ME
    refLength = -1;
    refHeader = ME.message;
end
end

function [seq, startPos, endPos, err] = getFastaSequenceAt(genomeBrowserDataMap, refSeqKey, startPos, endPos)
err='';
try
    obj= genomeBrowserDataMap.values({refSeqKey});
    obj = obj{:};
    
    startPos = startPos + 1;
    endPos = endPos + 1;
    if startPos < 1
        startPos = 1;
    end
    if isstruct(obj)
       totalLength = numel(obj.Sequence);
    else
        totalLength = obj.Length;
    end

    if endPos > totalLength
            endPos = totalLength;
    end 
        
    if isstruct(obj)
       seq = obj.Sequence(startPos:endPos); 
    else
        seq = obj.getSubSequence([startPos endPos]);
    end
catch ME
   err = ME.message;
   seq='';
   startPos =0;
   endPos = 0;
end    
end

function genomeBrowserDataMap = clearBioDataObject(genomeBrowserDataMap, bioDataKey)
if isKey(genomeBrowserDataMap, bioDataKey)
    genomeBrowserDataMap.remove(bioDataKey);
end
end

function [varName, varType] = getBioDataVariables()
vars = evalin('base', 'whos');

nameCell = {vars.name};
typeCell = {vars.class};

% % varFielName = fieldnames(vars);
% % % Find unwanted field names and remove from the vars struct
% % idx = findmatch(varFielName, {'name', 'class'});
% % vars = rmfield(vars, varFielName(~idx));
bdIdx = findmatch(typeCell,{'BioMap'});
varName = nameCell(bdIdx);
varType = typeCell(bdIdx);
end

function idx = findmatch(strs, lib)
% strs- cell array of strings
% lib - cell arrays of strings
% return indices of lib strings in strs
numLib = numel(lib);
idx = [];
for i = 1:numLib
     idx1 = strcmp(strs,lib{i});

     if isempty(idx)
         idx = idx1;
     else
        idx = idx | idx1;
     end
end
    
end

function bioMapKey = getGenomeDataMapKey(genomeBrowserDataMap, bioMapKey)
currentKeys = keys(genomeBrowserDataMap);
count = 1;
baseBioMapKey = bioMapKey;
while any(strcmp(currentKeys, bioMapKey))
    count = count+1;
    bioMapKey = sprintf('%s(%d)', baseBioMapKey, count);
end
end

function seqNames = getSequenceDictionary(mapKey)
seqNames = {};

if ~isempty(strfind(lower(mapKey), '.sam'))
    % sam file
    info = saminfo(mapKey);
else
    % bam file
    info = baminfo(mapKey);
end

if isfield(info, 'SequenceDictionary') && isfield(info.SequenceDictionary, 'SequenceName')
    seqNames = {info.SequenceDictionary.SequenceName};
end

end

function [seqNames, seqCount, errorId] = getScannedDictionary(mapKey)
seqNames = '';
seqCount = '';
errorId = '';

try
    if ~isempty(strfind(lower(mapKey), '.sam'))
        % sam file
        info = saminfo(mapKey, 'SCANDICTIONARY', true);
    else
        % bam file
        info = baminfo(mapKey, 'SCANDICTIONARY', true);
    end
    seqNames = {info.ScannedDictionary};
    seqCount = {info.ScannedDictionaryCount};
catch me
    errorId = me.identifier;
end

end




