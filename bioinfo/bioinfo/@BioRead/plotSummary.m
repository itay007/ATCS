function [hax, qs] = plotSummary(br,varargin)
% PLOTSUMMARY Summary plots of a BioRead object
%  PLOTSUMMARY(BR) generates a figure containing six plots that
%  present summary statistics of the data stored in the Illumina 1.8+ file
%  referenced by the BIOREAD object BR. 
%
% PLOTSUMMARY(...,'Encoding',ENC,...) determines how the characters
% encoding the scores in the FASTQ-file associated with BR should be
% interpreted. ENC can have the values 'Sanger', 'Illumina13', 
% 'Illumina15', 'Illumina18', and 'Solexa'. The default is 'Illumina18'.
%
% PLOTSUMMARY(...,'FilterLength', LEN,...) causes only the first
% LEN characters of each read to be used. LEN can be an integer greater
% than 0 or the empty array. In the latter case, no filtering is applied.
% The default is the empty array.
%
% PLOTSUMMARY(...,'QualityScoreThreshold', Q,...) results in
% reads with an average quality less than Q to be ignored. Q can have
% any scalar value. The default is -Inf, which causes all reads to be
% considered. 
%
% Note that if 'FilterLength' is set to LEN and 'QualityScoreThreshold' is
% set to Q, then a read is discarded if the average quality of the first
% LEN scores is less than Q.
%
% [H QC] = PLOTSUMMARY(...) returns a column vector H of handles to
% the axes in the generated figure and a QualityStatistics object QC that
% stores the data represented by the graphs.
%
% Example:
%
%  br = BioRead('SRR005164_1_50.fastq');
%  plotSummary(br,'Encoding','Illumina18','FilterLength',40,...
%                 'QualityScoreThreshold',5);


checkScalarInput(br);
details = getAdapterDetails(br);


if details.InMemory
    dataSource.Quality = {br.Quality{:}};
    dataSource.Sequence= {br.Sequence{:}};
    subset = NaN;
else
    dataSource = details.FileName;
    if details.IsSubset
        subset = details.SubsetIndex;
    else
        subset = NaN;
    end
end

qs = BioReadQualityStatistics(dataSource,'Subset',subset,varargin{:});
haxt = qs.plotSummary();
if nargout>0
    hax = haxt;
end


