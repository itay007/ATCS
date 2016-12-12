function write(obj,filename,varargin)
%WRITE writes the contents of a BioRead object.
%
%   WRITE(OBJ,FILENAME) writes the contents of a BioRead object to file
%   FILENAME. FILENAME is a string containing the base name of the file and
%   it can be prepended by an absolute or relative path. If the path is
%   missing the file is written in the same directory as the source file is
%   (when the object is indexed) or in the current directory (when the data
%   is in memory).
%
%   WRITE(...,'FORMAT',F) specifies the type of file format. Available
%   options are 'fasta','fastq'. Default is 'fasta' when the 'Quality'
%   property is empty, otherwise, default is 'fastq'.
%
%   WRITE(...,'OVERWRITE',TRUE) allows overwrite an existing file as long
%   as system file permissions are available. Default is FALSE. WRITE also
%   deletes the respective index file (.IDX) that becomes stale.

checkScalarInput(obj);
details = getAdapterDetails(obj);

if details.InMemory && isempty(obj.Quality)
    format = 'fasta';
else
    format = 'fastq';
end

% Parse optional PVPs and/or set defaults
[format, overwrite] = parse_inputs(format,varargin{:});

% Validate filename (if no file extension given one is added)
[outfilePath, outfileName, outfileExt] = fileparts(filename);
if ~isempty(outfileExt)
   if ~strcmpi(['.' format],outfileExt)
       error(message('bioinfo:BioRead:write:InvalidFileExtension'))
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
fileToCheck = [outfilePathName,'.',format,'.idx'];

if overwrite
    if strcmpi(filename,details.FileName)
        error(message('bioinfo:BioRead:write:OutputFileInvalid'))
    end    
    if exist(fileToCheck,'file')
        delete(fileToCheck)
        if exist(fileToCheck,'file')
            error(message('bioinfo:BioRead:write:CannotDeleteFile',fileToCheck))
        end
    end
else
    if exist(filename,'file')
        error(message('bioinfo:BioRead:write:OutputFileExists',filename))
    end
    if exist(fileToCheck,'file')
            warning(message('bioinfo:BioRead:write:StaleFile',fileToCheck))
    end  
end

if strcmp(details.FileFormat,'sam') || strcmp(details.FileFormat,'bam')
    bioinfoprivate.bamaccessmex('bam2bam',details.FileName,details.FileFormat,filename,format,uint32([details.SubsetIndex]))
    return;
end

if details.InMemory 
    if isempty(obj.Sequence) || (strcmp(format,'fastq') && isempty(obj.Quality))
        error(message('bioinfo:BioRead:write:EmptyProperties'))
    end
end    

fid = fopen(filename,'wt');
if fid<0
   error(message('bioinfo:BioRead:write:CannotOpenTextFile',filename))
end

blockSize = 1000;
try
    if strcmp(format,'fastq')
        for i = 1:blockSize:obj.NSeqs
            t = (get(getSubset(obj,i:min(i+blockSize,obj.NSeqs)),{'Header','Sequence','Quality'}));
            tmp = [t{1} t{2} t{3}]';
            fprintf(fid,'@%s\n%s\n+\n%s\n',tmp{:});
        end
    else
        for i = 1:blockSize:obj.NSeqs
            t = (get(getSubset(obj,i:min(i+blockSize,obj.NSeqs)),{'Header','Sequence'}));
            tmp = [t{1} t{2}]';
            fprintf(fid,'>%s\n%s\n',tmp{:});
        end
    end
catch ME
    fclose(fid);
    error(message('bioinfo:BioRead:write:CannotWriteToTextFile',filename))
end
fclose(fid);
end

function [format, overwrite] = parse_inputs(format,varargin)
% Parse input PV pairs.

% default values
overwrite = false;

if rem(nargin, 2) == 0
    error(message('bioinfo:BioRead:write:IncorrectNumberOfArguments', mfilename));
end
okargs = {'format', 'overwrite'};
for j=1:2:nargin-1
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1  % format
            [~,format] = bioinfoprivate.optPartialMatch(pval,{'fasta','fastq'}, okargs{k}, mfilename); 
         case 2 % overwrite
            overwrite = bioinfoprivate.opttf(pval,okargs{k},mfilename);
    end
end
end