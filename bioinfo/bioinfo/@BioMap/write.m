function write(obj,filename,varargin)
%WRITE writes the contents of a BioMap object.
%
%   WRITE(OBJ,FILENAME) writes the contents of a BioMap object to file
%   FILENAME. FILENAME is a string containing the base name of the file and
%   it can be prepended by an absolute or relative path. If the path is
%   missing the file is written in the same directory as the source file is
%   (when the object is indexed) or in the current directory (when the data
%   is in memory).
%
%   WRITE(...,'FORMAT',F) specifies the type of file format. Available
%   options are 'fasta','fastq','bam','sam'. Default is 'bam'.
%
%   WRITE(...,'OVERWRITE',TRUE) allows overwrite an existing file as long
%   as system file permissions are available. Default is FALSE. WRITE also
%   deletes index files (.IDX,.BAI,.LINEARINDEX) and ordered files
%   (.ORDERED.BAM, .ORDERED.SAM) that become stale.

%  Copyright 2012 The MathWorks, Inc.
	
checkScalarInput(obj);
details = getAdapterDetails(obj);

% Parse optional PVPs and/or set defaults
[format, overwrite] = parse_inputs(varargin{:});

% Validate filename (if no file extension given one is added)
[outfilePath, outfileName, outfileExt] = fileparts(filename);
if ~isempty(outfileExt)
   if ~strcmpi(['.' format],outfileExt)
       error(message('bioinfo:BioMap:write:InvalidFileExtension'))
   end
end
outfileExt = ['.' format];

% If no file path specified, we use the same location where the source is,
% for objects in memory we use the current path
if isempty(outfilePath) && ~isempty(details.FileName)
    outfilePath = fileparts(details.FileName);
end

outfilePathName = fullfile(outfilePath,outfileName);
filename = [outfilePathName, outfileExt];

switch format
    case 'sam',
        filesToCheck = {[outfilePathName,'.sam.idx']
                        [outfilePathName,'.ordered.sam']
                        [outfilePathName,'.ordered.sam.idx']};
    case 'bam',
        filesToCheck = {[outfilePathName,'.bam.bai']
                        [outfilePathName,'.bam.linearindex']
                        [outfilePathName,'.ordered.bam']
                        [outfilePathName,'.ordered.bam.bai']
                        [outfilePathName,'.ordered.bam.linearindex']};
    case 'fasta',
        filesToCheck = {[outfilePathName,'.fasta.idx']};
    case 'fastq',
        filesToCheck = {[outfilePathName,'.fastq.idx']};
end
            
if overwrite
    if strcmpi(filename,details.FileName)
        error(message('bioinfo:BioMap:write:OutputFileInvalid'))
    end
    for i=1:numel(filesToCheck)
        if exist(filesToCheck{i},'file')
            delete(filesToCheck{i})
            if exist(filesToCheck{i},'file')
                error(message('bioinfo:BioMap:write:CannotDeleteFile',filesToCheck{i}))
            end
        end
    end
else
    if exist(filename,'file')
        error(message('bioinfo:BioMap:write:OutputFileExists',filename))
    end
    for i=1:numel(filesToCheck)
        if exist(filesToCheck{i},'file')
            warning(message('bioinfo:BioMap:write:StaleFile',filesToCheck{i}))
        end
    end    
end
if details.InMemory
    if isempty(obj.SequenceDictionary) || isempty(obj.Header) || ...
       isempty(obj.Flag) || isempty(obj.Reference) || isempty(obj.Start) || ...
       isempty(obj.MappingQuality) || isempty(obj.Signature) || ...
       isempty(obj.Sequence) || isempty(obj.Quality)
        error(message('bioinfo:BioMap:write:EmptyProperties'))
    end
    if strcmp(format,'sam')
        tname = filename;
    else
        tname = [tempname '.sam'];
    end
    fid = fopen(tname,'wt');
    if fid<0
       error(message('bioinfo:BioMap:write:CannotOpenTextFile',tname))
    end
    try
       for i = 1:numel(obj.SequenceDictionary)
           fprintf(fid,'@SQ\tSN:%s\tLN:536870911\n',obj.SequenceDictionary{i});
           % SN is required, however BioMap does not have this information
           % so we give the largest possible number in the format
           % specification
       end
       for i = 1:obj.NSeqs
           fprintf(fid,'%s\t%d\t%s\t%d\t%d\t%s\t*\t0\t0\t%s\t%s\n',obj.Header{i},obj.Flag(i),obj.Reference{i},obj.Start(i),obj.MappingQuality(i),obj.Signature{i},obj.Sequence{i},obj.Quality{i});
       end
    catch ME
       fclose(fid);
       error(message('bioinfo:BioMap:write:CannotWriteToTextFile',tname)) 
    end
    fclose(fid);
    if ~strcmp(format,'sam')
        bioinfoprivate.bamaccessmex('bam2bam',tname,'sam',filename,format,uint32([])) 
        delete(tname)
    end
else
    bioinfoprivate.bamaccessmex('bam2bam',details.FileName,details.FileFormat,filename,format,uint32([details.SubsetIndex]))
end

end

function [format, overwrite] = parse_inputs(varargin)
% Parse input PV pairs.

% default values
format = 'bam';
overwrite = false;

if rem(nargin, 2) ~= 0
    error(message('bioinfo:BioMap:write:IncorrectNumberOfArguments', mfilename));
end
okargs = {'format', 'overwrite'};
for j=1:2:nargin-1
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1  % format
            [~,format] = bioinfoprivate.optPartialMatch(pval,{'sam','bam','fasta','fastq'}, okargs{k}, mfilename); 
         case 2 % overwrite
            overwrite = bioinfoprivate.opttf(pval,okargs{k},mfilename);
    end
end
end