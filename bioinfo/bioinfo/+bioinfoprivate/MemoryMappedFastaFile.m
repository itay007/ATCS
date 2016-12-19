classdef  MemoryMappedFastaFile
%MEMORYMAPPEDFASTAFILE Class that allows indexed access to a FASTA file.
%
%  A MEMORYMAPPEDFASTAFILE object stores information that allows the
%  efficient random read access to a FASTA formatted file. The FASTA
%  formatted file must contain a single entry and equally sized rows.
%
%  MemoryMappedFastaFile properties:
%    FileName    - Location and name of the FASTA formatted file.
%    Header      - Header in the FASTA formatted file.
%    Length      - Length of the sequence (including valid positions
%                  not in the FASTA file). 
%    Range       - Valid positions in the object.
%    MappedRange - Positions available from the FASTA file.        
%    Offset      - Offset between the positions in the object and the
%                  positions in the FASTA file. 
%
%  MemoryMappedFastaFile methods:
%    MemoryMappedFastaFile - creates a MemoryMappedFastaFile object.
%    getSubSequence        - returns a subsequence from the object.
%    Sequence              - subscripted indexing into the object.
%
%  Example:
%    ref = MemoryMappedFastaFile('Homo_sapiens.GRCh37.56.dna.chromosome.1.fa')
%    seq = getSubSequence(ref,[10000000,10100000]);
%    seq = ref.Sequence(10000000:10100000);
%
%  See also fastaread, BioIndexedFile, MemoryMappedFastaFile.getSubSequence,
%          MemoryMappedFastaFile.Sequence.

%   Copyright 2010-2012 The MathWorks, Inc.


    properties (SetAccess = private, GetAccess = public)
%FILENAME  Location and name of the FASTA formatted file.
        FileName
    end
    properties (SetAccess = public, GetAccess = public)
%HEADER Header in the FASTA formatted file.        
        Header
    end
    properties (Dependent = true, SetAccess = private, GetAccess = public)
%LENGTH Length of the sequence (including valid positions not in the FASTA file).        
        Length
    end     
    properties (SetAccess = private, GetAccess = public) 
%RANGE Valid positions in the object.        
        Range
    end
    properties (Dependent = true, SetAccess = private, GetAccess = public)
%MAPPEDRANGE Positions available from the FASTA file.        
        MappedRange
    end  
    properties (SetAccess = private, GetAccess = public)
%OFFSET Offset between the positions in the object and the positions in the file.        
        Offset
    end
    
    properties (SetAccess = private, GetAccess = private)
        OffsetToSequence
        RowLengthInFile
        RowLengthInSequence
        NucleotideFileLength
    end    

    methods
        
        function l = get.Length(obj)
            l = diff(obj.Range)+1;
        end
        
        function r = get.MappedRange(obj)
            r = [max(obj.Range(1),obj.Offset+1) min(obj.Range(2),obj.Offset+obj.NucleotideFileLength)];
        end
        
        function obj = MemoryMappedFastaFile(fastafilename,offset,range)
%MEMORYMAPPEDFASTAFILE creates a MemoryMappedFastaFile object.
%
%  MEMORYMAPPEDFASTAFILE(FASTAFILE) constructs a MemoryMappedFastaFile
%  object that indexes the contents of a FASTA formatted file. The FASTA
%  formatted file must contain a single entry and equally sized rows.            
%
%  MEMORYMAPPEDFASTAFILE(FASTAFILE, OFFSET) indicates the offset of the
%  sequence in the FASTA formatted file relative with the actual reference.
%  OFFSET is an integer. OFFSET defaults to 0, which indicates that
%  positions in the FASTA formatted file are the same as the real positions
%  in the reference. Valid positions in the object ('Range' property) is
%  set to [1 L+OFFSET]. Where L is the length of the FASTA formatted file
%  (in nucleotides). Queries to valid positions but not available in the
%  FASTA formatted file return an 'N' symbol.
%
%  MEMORYMAPPEDFASTAFILE(FASTAFILE, OFFSET, RANGE) Indicates the valid
%  positions in the object. RANGE is a 1-by-2 vector with positive
%  integers. Queries to valid positions but not available in the FASTA
%  formatted file return an 'N' symbol.
%
%  Example:
%    ref = MemoryMappedFastaFile('Homo_sapiens.GRCh37.56.dna.chromosome.1.fa')
%    seq = getSubSequence(ref,[10000000,10100000]);
%    seq = ref.Sequence(10000000:10100000);
%
%  See also fastaread, BioIndexedFile, MemoryMappedFastaFile.getSubSequence
            lookAheadLines = 1000; 
            lookBackLines = 10;
            if isempty(fileparts(fastafilename)) && exist(fastafilename,'file')
              fastafilename = which(fastafilename);
            end
            obj.FileName = fastafilename;
            fid = fopen(obj.FileName);
            if fid<0
                if exist(fastafilename,'file')
                   error(message('bioinfo:MemoryMappedFastaFile:MemoryMappedFastaFile:CannotOpenFile'))
                else
                   error(message('bioinfo:MemoryMappedFastaFile:MemoryMappedFastaFile:CannotFindFile'))
                end
            else
                c = onCleanup(@()fclose(fid));
            end
            % Analyze first line, must be header:
            fileLine = fgetl(fid);
            if isempty(fileLine) || fileLine(1)~='>'
                error(message('bioinfo:MemoryMappedFastaFile:MemoryMappedFastaFile:InvalidHeader'))
            end
            obj.Header =  regexprep(fileLine,'^>','');
            % Anlayze file position and sequence position for the next
            % lookAheadLines lines or EOF
            count = 1;
            filePos = zeros(lookAheadLines,1);
            filePos(1) = ftell(fid);
            seqPos = zeros(lookAheadLines,1);
            while count<lookAheadLines 
                fileLine = fgetl(fid);
                if isempty(fileLine) || fileLine(1)==-1
                    break
                end
                if fileLine(1)=='>'
                    error(message('bioinfo:MemoryMappedFastaFile:MemoryMappedFastaFile:MultipleEntries'))
                end
                count = count + 1;
                filePos(count) = ftell(fid);
                seqPos(count) = numel(regexprep(fileLine,'\W',''));
            end
            if count == 2
                fileOffsets = seqPos(2); 
                seqOffsets = seqPos(2);
            elseif isempty(fileLine) || fileLine(1)==-1
                fileOffsets = unique(diff(filePos(1:count-1)));
                seqOffsets = unique(seqPos(2:count-1));
            else
                fileOffsets = unique(diff(filePos));
                seqOffsets = unique(seqPos(2:end));
            end
            if numel(fileOffsets)>1 || numel(seqOffsets)>1
                error(message('bioinfo:MemoryMappedFastaFile:MemoryMappedFastaFile:UnequalSizedRows'))
            end
            obj.OffsetToSequence = filePos(1);
            obj.RowLengthInFile = fileOffsets;
            obj.RowLengthInSequence = seqOffsets;
            % Analyze last line and figure out the Sequence Length
            status = fseek(fid,0,'eof');
            if status~= 0 
                error(message('bioinfo:MemoryMappedFastaFile:MemoryMappedFastaFile:CannotReadLastLine'))
            end
            eof = ftell(fid);
            status = fseek(fid,-obj.RowLengthInFile*lookBackLines,'eof');
            if status~= 0 
                status = fseek(fid,1,'bof');
                if status~= 0
                    error(message('bioinfo:MemoryMappedFastaFile:MemoryMappedFastaFile:CannotReadLastLine'))
                end
            end            
            fgetl(fid);
            filePos = ftell(fid);
            nCharsLastLine = 0;
            while filePos~=eof
               filePos = ftell(fid);
               fileLine = fgetl(fid); 
               if isempty(fileLine) || fileLine(1)==-1
                   break
               end
               nChars = numel(regexprep(fileLine,'\W',''));
               if nChars>0 
                   nCharsLastLine = nChars;
                   filePosLastLine = filePos;
                   if  nCharsLastLine~=obj.RowLengthInSequence
                       break
                   end
               end
            end
            if nCharsLastLine==0
                 error(message('bioinfo:MemoryMappedFastaFile:MemoryMappedFastaFile:CannotReadLastLine'))
            end
            numLines = (filePosLastLine-obj.OffsetToSequence)/obj.RowLengthInFile;
            if  rem(numLines,1)
                error(message('bioinfo:MemoryMappedFastaFile:MemoryMappedFastaFile:UnequalSizedRows'))
            end
            obj.NucleotideFileLength = numLines .* obj.RowLengthInSequence + nCharsLastLine;
            
            if nargin == 1
                obj.Offset = 0;
                obj.Range = [1 obj.NucleotideFileLength];
            else % nargin>1
                if ~isnumeric(offset) || ~isscalar(offset) || rem(offset,1)
                     error(message('bioinfo:MemoryMappedFastaFile:MemoryMappedFastaFile:InvalidOffset'))
                end
                obj.Offset = offset;
                obj.Range = [1 obj.NucleotideFileLength] + obj.Offset;
                if nargin==2
                    if obj.Range(2)<1
                        error(message('bioinfo:MemoryMappedFastaFile:MemoryMappedFastaFile:CorruptingOffset'))
                    end
                    obj.Range(1) = 1;
                else % nargin==3
                    if ~isnumeric(range) || numel(range)~=2 || any(rem(range,1)) || range(1)>range(2) || range(1)<1
                       error(message('bioinfo:MemoryMappedFastaFile:MemoryMappedFastaFile:InvalidRange'))
                    end
                    if range(1)>obj.Range(2) || range(2)<obj.Range(1)
                       error(message('bioinfo:MemoryMappedFastaFile:MemoryMappedFastaFile:InvalidRangeOffset'))
                    end
                    obj.Range = range(:)';
                end
            end

        end
        
        function seq = getSubSequence(obj,X)
        %GETSUBSEQUENCE Returns a subsequence from the object.
        %
        %  GETSUBSEQUENCE(OBJ,X) Returns a range of the sequence from the
        %  memory mapped FASTA file starting at sequence position X(1) and
        %  ending at sequence position X(2). X is a two element row vector
        %  with positive integers such that X(1)<= X(2).
        %
        %  Example:
        %    ref = MemoryMappedFastaFile('Homo_sapiens.GRCh37.56.dna.chromosome.1.fa')
        %    seq = getSubSequence(ref,[10000000,10100000]);
        %
        %  See also  MemoryMappedFastaFile, MemoryMappedFastaFile.Sequence
            
            if  ~ismatrix(X) || size(X,2)~=2 || size(X,1)>1 || any(X(:)<1) || ...
                any(rem(X(:),1)~=0) || any(diff(X,[],2)<0) 
                error(message('bioinfo:MemoryMappedFastaFile:getSubSequence:InvalidRange'));
            end
            
            if X(1)<obj.Range(1) || X(2)>obj.Range(2)
                error(message('bioinfo:MemoryMappedFastaFile:getSubSequence:RangeOutOfBounds'));
            end

            fid = fopen(obj.FileName);
            if fid<0
                error(message('bioinfo:MemoryMappedFastaFile:getSubSequence:CannotOpenFile'))
            else
                c = onCleanup(@()fclose(fid));
            end       
            seq = readSequenceFromObject(obj,X,fid);
        end
        
        function seq = Sequence(obj,i)
        %SEQUENCE Subscripted indexing into a memory mapped FASTA file.
        %
        %  SEQUENCE(OBJ,I) One dimensional subscripted indexing into a
        %  memory mapped FASTA file. I is a numeric vector with positive
        %  integers within the range of valid positions. There is a
        %  one-to-one relationship between the values in I and the output,
        %  order and quantity are preserved, despite repeated values in I.
        %  However, access to the FASTA file is as efficient as possible by
        %  minimizing the number of file seeks and read operations.
        %
        %  Example:
        %    ref = MemoryMappedFastaFile('Homo_sapiens.GRCh37.56.dna.chromosome.1.fa')
        %    seq = ref.Sequence(10000000:10100000);
        %
        %  See also  MemoryMappedFastaFile, MemoryMappedFastaFile.getSubSequence
        
            maxSkip = 100;
            if ~isnumeric(i) || any(i<1) || any(rem(i,1))
                 error(message('bioinfo:MemoryMappedFastaFile:Sequence:nonNumeric'))
            end
            if isscalar(i)
                seq = getSubSequence(obj,[i,i]);
            elseif isempty(i)
                seq = '';
            elseif all(diff(i)==1)
                seq = getSubSequence(obj,[i(1),i(end)]);
            else
                fid = fopen(obj.FileName);
                if fid<0
                    error(message('bioinfo:MemoryMappedFastaFile:Sequence:CannotOpenFile'))
                else
                    c = onCleanup(@()fclose(fid));
                end
                seq(1,numel(i)) = ' ';
                [u,~,i] = unique(i(:));
                blocks = [find([1;diff(u)>maxSkip]);numel(u)+1];
                for idx = 1:(numel(blocks)-1)
                    r = readSequenceFromObject(obj,[u(blocks(idx)),u(blocks(idx+1)-1)],fid);
                    li = i>=blocks(idx) & i<=blocks(idx+1)-1;
                    seq(li) = r(u(i(li))-u(blocks(idx))+1);
                end
            end            
        end
    end
    
    methods (Access = private)
        
        function seq = readSequenceFromObject(obj,X,fid)
            i = max(obj.MappedRange(1),X(1));
            j = min(obj.MappedRange(2),X(2));
            if j>=i
                seq = [ char(blanks(i-X(1))+'.') upper(readSequenceFromFastaFile(obj,i-obj.Offset,j-obj.Offset,fid))  char(blanks(X(2)-j)+'.') ];
            else
                seq = [ char(blanks(i-X(1))+'.') char(blanks(X(2)-j)+'.') ];
            end
        end
            
        function seq = readSequenceFromFastaFile(obj,i,j,fid)
            % i is in fileLine ib:
            ib = floor(double(i-1)/obj.RowLengthInSequence)+1;
            % start of line ib in the file is:
            fpos_ib = (ib-1)*obj.RowLengthInFile+obj.OffsetToSequence;
            % start of line ib in the sequence is:
            spos_ib = (ib-1)*obj.RowLengthInSequence+1;
            % extra chars in the first line
            offset_first_line = i - spos_ib;
            % length of query
            le = j-i+1;
            % number of complete lines expected to read:
            nl = floor(double(offset_first_line+le-1)/obj.RowLengthInSequence);
            % number of bytes to read
            nb = (obj.RowLengthInFile-obj.RowLengthInSequence)*nl + le;
            % file offset
            fo = fpos_ib + offset_first_line;
            % move to file position
            status = fseek(fid,fo,'bof');
            if status~= 0 
                error(message('bioinfo:MemoryMappedFastaFile:MemoryMappedFastaFile:CannotMoveToFilePosition'))
            end
            [seq count] = fscanf(fid,'%c',nb);
            if count ~= nb
                   error(message('bioinfo:MemoryMappedFastaFile:MemoryMappedFastaFile:CannotReadRange'))
            end
            seq(seq<'A')=[];  
            seq(seq>'z')=[]; 
        end
    end
end
