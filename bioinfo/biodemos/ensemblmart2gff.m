function out = ensemblmart2gff(inp,out)
%ENSEMBLMART2GFF converts from ENSEMBL's BioMart tables to GFF
%
% ENSEMBLMART2GFF is a helper function for MBDSEQDEMO and RNASEQDEDEMO.
%
% ENSEMBLMART2GFF(INPUTFILENAME,OUTPUTFILENAME) converts a TSV table
% obtained with the ENSEMBL's BioMart service into a GFF formatted file. 
% Data fields are translated as follows:
%
%   GFF Fields      ENSEMBL Attributes
%   seqname         chromosome_name or 'unknown' if attribute is missing
%   source          gene_biotype or 'EnsEMBL' if attribute is missing
%   feature         associated_gene_name or 'gene' if attribute is missing
%   start           gene_start_(bp)
%   end             gene_end_(bp)
%   score           '0.0'
%   strand          Strand or '+' if attribute is missing
%   frame           '.'
%   attribute        all other existing ENSEMBL Attributes
% 
% OUTFILENAME = ENSEMBLMART2GFF(INPUTFILENAME) figures out automatically
% the name of the GFF file by using the same path and name as INPUTFILENAME
% but changing the extension to GFF.
%
% Note: This utility function is limited for simple gene lists. For
% structural information you can directly download GFF and GTF formatted
% files from the ENSEMBL's ftp site or using ENSEMBL's BioMart service.

%   Copyright 2012 The MathWorks, Inc.

 fid = fopen(inp);

 if fid<0
    error('bioinfo:ensemblmart2gff:invalidFile','Cannot open input file.')
 end
     
 c = onCleanup(@()fclose(fid));
 h = fgetl(fid);
 a = textscan(h,'%s','delimiter','\t');
 a = a{1};
 n = numel(a);
 
 GFFseqname = find(strncmpi(a,'Chromosome Name',numel('Chromosome Name')));
 if numel(GFFseqname)~=1
     GFFseqname = 'unknown';
 end
 
 GFFsource = find(strncmpi(a,'Gene Biotype',numel('Gene Biotype')));
 if numel(GFFsource)~=1
     GFFsource = 'EnsEMBL';
 end
 
 GFFfeature = find(strncmpi(a,'Associated Gene Name',numel('Associated Gene Name')));
 if numel(GFFfeature)~=1
     GFFfeature = 'gene';
 end
 
 GFFstart = find(strncmpi(a,'Gene Start (bp)',numel('Gene Start (bp)')));
 if numel(GFFstart)~=1
     error('bioinfo:ensemblmart2gff:invalidStart','Column ''Gene Start (bp)'' must be available in the input file.')
 end
 
 GFFend = find(strncmpi(a,'Gene End (bp)',numel('Gene End (bp)')));
 if numel(GFFend)~=1
     error('bioinfo:ensemblmart2gff:invalidEnd','Column ''Gene End (bp)'' must be available in the input file.')
 end
 
 GFFscore = '0.0';
 
 GFFstrand = find(strncmpi(a,'Strand',numel('Strand')));
 if numel(GFFstrand)~=1
     GFFstrand = '+';
 end
 
 GFFframe = '.';
 
 GFFattribute = true(n,1);
 if isnumeric(GFFseqname) 
     GFFattribute(GFFseqname) = false; 
 end
 if isnumeric(GFFsource)  
     GFFattribute(GFFsource) = false;  
 end
 if isnumeric(GFFfeature) 
     GFFattribute(GFFfeature) = false; 
 end
 GFFattribute(GFFstart) = false;
 GFFattribute(GFFend) = false;
 if isnumeric(GFFstrand) 
     GFFattribute(GFFstrand) = false; 
 end
 GFFattribute = find(GFFattribute);
 
 att = regexprep(a,'\s','_');
 att = strcat(att(GFFattribute),' "%s"; ');
 att = sprintf('%s ',att{:});

 if nargin==1 
    [path,name,ext] = fileparts(inp);
    if any(strncmpi({'.GFF','.GTF'},ext,4))
         error('bioinfo:ensemblmart2gff:invalidExtension','Input file has invalid extension.')
    end
    out = fullfile(path,[name '.gff']);
 end
    
 fid2 = fopen(out,'wt');
 c2 = onCleanup(@()fclose(fid2));
  
 while true
     d = textscan(fid,'%s',n,'Delimiter','\t');
     d = d{1};
     if isempty(d)
         break
     end
     strand = strrep(regexprep(sel(d,GFFstrand),'-1','-'),'1','+');
     fprintf(fid2,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
             sel(d,GFFseqname),sel(d,GFFsource),sel(d,GFFfeature),...
             d{GFFstart},d{GFFend},sel(d,GFFscore),strand,...
             sel(d,GFFframe),sprintf(att,d{GFFattribute}));
 end
 
    function str = sel(d,i)
        if isnumeric(i)
            str = d{i};
        else
            str = i;
        end
    end

end
 
 