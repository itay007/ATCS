classdef BioReadQualityStatistics
	%BIOREADQUALITYSTATISTICS quality statistics from short-read sequences.
    %
    %  The BIOREADQUALITYSTATISTICS class contains quality statistics data
    %  from short-read sequences and provides a standard set of quality
    %  control plots for such data.   
    %
    %  Construct a BioReadQualityStatistics object from short-read sequence
    %  data stored in FASTQ, SAM, or BAM files. Perform data quality
    %  analyses using the object's methods to generate several quality
    %  control plots regarding average quality score for each base
    %  position, average quality score distribution, read count percentage
    %  for each base position, percentage of G and C nucleotides for each
    %  base position, G and C content distribution, and all nucleotide
    %  distribution. The object lets parse a sequence file without creating
    %  a BioRead object and interact with the quality data in order to
    %  compare different data sets or filtering options and create
    %  customized plots.
	%
	%  BioReadQualityStatistics methods:
    %  plotPerPositionCountByQuality - plot percentages of reads with scores.
    %  plotPerPositionGC       - plot percentages of G or C nucleotides.
    %  plotPerPositionQuality  - plot Phred score distributions.
    %  plotPerSequenceGC       - plot G or C nucleotide distribution.
    %  plotPerSequenceQuality  - plot distribution of average quality scores.
    %  plotSummary             - plot all summary statistics in one figure.
    %  plotTotalGC             - plot distribution of all nucleotides.
    %
    % Example:
    %      qc = BioReadQualityStatistics('SRR005164_1_50.fastq',...
    %              'Encoding','Illumina18','FilterLength',40,...
    %              'QualityScoreThreshold',5); 
    %      plotSummary(qc)
    %
	%  See also BIOMAP, BIOREAD.
	
	%  Copyright 2012-2013 The MathWorks, Inc.    
    
    properties (GetAccess = public, SetAccess = private)
        FileName;
        FileType;
        Encoding;
        CharOffset;
        NumberOfReads;
        MaxReadLength;
        MinEncodingPhred;
        MaxEncodingPhred;
        SkipPhred;
        PerSeqAverageQualityDist;
        PerPosQualities;
        PerSeqGCDist;
        PerPosBaseDist;
        Name;
        
        MaxScore;
        MinScore;
        
        FilterLength;
        QualityScoreThreshold;
        
        Subset;
        
    end
    
    properties (Constant, Access = private)
        NoSkip = 256;
        NoFilterLength = -1;
        NoSubset = -1;
        ValidEncodings = {'Sanger','Illumina13','Illumina15' ,'Illumina18','Solexa'};
        ValidFileTypes = {'FASTQ', 'SAM','BAM', 'Memory'};
    end
    
    methods (Access = public)
        function obj = BioReadQualityStatistics(varargin)
            %BIOREADQUALITYSTATISTICS Create a QualityStatistics Object
            %
            % QS = BIOREADQUALITYSTATISTICS(FNAME) creates a QualityStatistics object from
            % the data stored in the Illumina 1.8+ file FNAME.
            %
            % QS = BIOREADQUALITYSTATISTICS(...,'Encoding',ENC,...) determines how the file
            % should be interpreted. ENC can have the value 'Sanger', 'Illumina13',
            % 'Illumina15', 'Illumina18', or 'Solexa'. The default is 'Illumina18'.
            %
            % QS = BIOREADQUALITYSTATISTICS(...,'FilterLength', LEN,...) causes only the first
            % LEN characters of each read to be used. LEN can be an integer greater
            % than 0 or the empty array. In the latter case, no filtering is applied.
            % The default is the empty array.
            %
            %  QS = BIOREADQUALITYSTATISTICS(...,'QualityScoreThreshold', Q,...) results in
            %  reads with an average quality less than Q to be ignored. Q can have
            %  any value. The default is -Inf, which causes all reads to be considered.
            %
            %  QS = BIOREADQUALITYSTATISTICS(...,'Subset', S,...) causes only reads in 
            %  the vector S to be used. For example, if S = [1,2,4], then only 
            %  the first, second and fourth reads will be read. S can be NaN, in 
            %  which case all reads are considered, an empty array, in which case
            %  no reads or considered, of a strictly increasing array of integers.
            %  The default is NaN.
            %
            % Note that if 'FilterLength' is set to LEN and 'QualityScoreThreshold' is
            % set to Q, then a read is discarded if the average quality of the first
            % LEN characters is less than Q.
            %
            % Example:
            %      qc = BioReadQualityStatistics('SRR005164_1_50.fastq',...
            %                             'Encoding','Illumina18','FilterLength',40,...
            %                             'QualityScoreThreshold',5);
            
            
            
            ip = inputParser();
            ip.StructExpand = false;
            ip.addRequired('fname');
            ip.addParamValue('Encoding','Illumina18');
            ip.addParamValue('FilterLength',[]);
            ip.addParamValue('QualityScoreThreshold',-Inf);
            ip.addParamValue('Name','');
            ip.addParamValue('FileType','');
            ip.addParamValue('Subset',NaN);
            
            ip.parse(varargin{:});
            
            obj.Name = ip.Results.Name;
            obj.Encoding = validatestring(ip.Results.Encoding,BioReadQualityStatistics.ValidEncodings);
            [obj.MinEncodingPhred, obj.MaxEncodingPhred,...
                obj.CharOffset, obj.SkipPhred] = BioReadQualityStatistics.getStandardScales(obj.Encoding);
            obj.FilterLength = ip.Results.FilterLength;
            obj.QualityScoreThreshold = ip.Results.QualityScoreThreshold;
            obj.Subset = ip.Results.Subset;
            
            assert(isValidSubset(obj.Subset),'bioinfo:BioReadQualityStatistics:InvalidSubset',...
                   'Subset Indices must be a strictly increasing set of integers that are greater than 0.'); 
            
            
            if ischar(ip.Results.fname)
                [obj.FileName, obj.FileType, deltaFromFile] = setUpInfoForFile(ip,obj.CharOffset);
                dataSource = obj.FileName;
            else
                obj.FileName = ''; obj.FileType = 'MEMORY'; deltaFromFile = obj.CharOffset;
                dataSource = ip.Results.fname;
            end
            
            
            cTypeNum = BioReadQualityStatistics.getTypeNumber(obj.FileType);
            
            cSkip = ...
                BioReadQualityStatistics.giveValueIfEmpty(obj.SkipPhred,BioReadQualityStatistics.NoSkip);
            cFilterLength = ...
                BioReadQualityStatistics.giveValueIfEmpty(obj.FilterLength,BioReadQualityStatistics.NoFilterLength);
            
                        
            if isnan(obj.Subset)
                cSubsetLen = BioReadQualityStatistics.NoSubset;
                cSubset    = 0; % Not Used
            else 
                cSubsetLen = numel(obj.Subset);
                cSubset = obj.Subset-1; %For C style indexing
            end
                
            
            
            [obj.PerSeqAverageQualityDist, ...
                obj.PerPosQualities,...
                obj.PerSeqGCDist,...
                obj.PerPosBaseDist,...
                obj.NumberOfReads,...
                obj.MaxReadLength] = bioinfoprivate.scanQualities(...
                int64(cTypeNum),...
                dataSource,...
                int64(obj.MinEncodingPhred),...
                int64(obj.MaxEncodingPhred),...
                int64(deltaFromFile),...
                int64(cSkip),...
                int64(cFilterLength),...
                double(obj.QualityScoreThreshold),...
                int64(cSubsetLen),...
                int64(cSubset));
            
            [obj.MinScore, obj.MaxScore] = obj.getMinMaxScore();
            
        end
        
        
        
        
        hax = plotPerPositionQuality(qc,h_ax);
        hax = plotPerSequenceQuality(qc,hax);
        hax = plotPerSequenceGC(qc,hax);
        hax = plotTotalGC(qc,hax);
        hax = plotPerPositionGC(qc,hax);
        hax = plotPerPositionCountByQuality(qc,hax);
        hf = plotSummary(data);
        
        
    end
    
    methods (Access = private)
        function [mn, mx] = getMinMaxScore(obj)
            mn = NaN; mx = NaN;
            if obj.NumberOfReads == 0
                return;
            end
            range = find(any(obj.PerPosQualities,2));
            if isempty(range)
                return;
            end
            mn = obj.MinEncodingPhred + (range(1)  -1);
            mx = obj.MinEncodingPhred + (range(end)-1);
        end
        
        
        function lim = phredAxisLim(obj)
            if isnan(obj.MinPhred)
                lim = 0:40;
                return;
            end
            lim = min(0,obj.MinPhred):max(40,obj.MaxPhred);
        end
        
        
    end
    methods (Static = true, Access = private)
        function [mn, mx, offset, skip] = getStandardScales(filetype)
            switch filetype
                case 'Sanger'
                    mn = 0; mx = 93; offset = 33; skip = [];
                case 'Illumina13'
                    mn = 0; mx = 62; offset = 64; skip = [];
                case 'Illumina15'
                    mn = 2; mx = 62; offset = 64;  skip = 2;
                case 'Illumina18'
                    mn = 0; mx = 62; offset = 33; skip = [];
                case 'Solexa'
                    mn = -5; mx = 62; offset = 64; skip = [];
                otherwise
                    error('bioinfo:BioReadQualityStatistics:InvalidEncoding','Not a valid encoding.');
            end
        end
        
        function val = giveValueIfEmpty(val,newval)
            if isempty(val)
                val = newval;
            end
        end
        
        function typeNum = getTypeNumber(str)
            switch str
                case 'FASTQ'
                    typeNum = 1;
                case 'SAM'
                    typeNum = 2;
                case 'BAM'
                    typeNum = 3;
                case 'MEMORY'
                    typeNum = 4;
                otherwise
                    error('bioinfo:BioReadQualityStatistics:InvalidFileType','Not a valid File type.');
            end
        end
    end
end



function [FileName, FileType, deltaFromFile] = setUpInfoForFile(ip,CharOffset)
fp = which(ip.Results.fname);
if isempty(fp)
    FileName = ip.Results.fname;
else
    FileName = fp;
end
[~,~,ext] = fileparts(FileName);
if ~isempty(ip.Results.FileType)
    FileType = ip.Results.FileType;
elseif ~isempty(ext)
    FileType = validatestring(upper(ext(2:end)),BioReadQualityStatistics.ValidFileTypes);
else
    error('bioinfo:BioReadQualityStatistics:UnknownFileType','Unable to determine file type.');
end
FileType = upper(FileType);
if isequal(FileType,'SAM') || isequal(FileType,'BAM')
    deltaFromFile = CharOffset-33;  %SAMTOOLS handles the shifts
else
    deltaFromFile = CharOffset;
end
end


function valid = isValidSubset(x)

if isempty(x) || (isscalar(x)&&isnan(x))
    valid = true;
    return;
end
valid = all((rem(x,1) == 0) &  (x>0)) && all(diff(x)>=1);
end

