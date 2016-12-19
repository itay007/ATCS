function ngsbrowser(bioMapVar)
%NGSBROWSER browse short-read sequence alignment.
%
%   NGSBROWSER opens an empty session of the browser. From the browser, you
%   can load:
%
%      * Alignment data from a SAM-formatted file or from a BioMap
%        in the MATLAB Workspace;
%      * Reference sequence from a FASTA file;
%      * Annotation features from a GFF or GTF file. 
%
%   Example:
%       b = BioMap('ex1.sam');
%       ngsbrowser
%  
%   %   From the browser, select File->Import Alignment Data from MATLAB Workspace
%   %   to load b into the browser.
% 
%   See also BIOMAP, SEQALIGNVIEWER, SEQVIEWER. 
 
%    Copyright 2010-2012 The MathWorks, Inc. 

%   Undocumented:
%   NGSBROWSER(BM) loads the short read sequence alignments stored in
%   the BioMap object BM.

if ~usejava('jvm') || ~usejava('awt') || ~usejava('swing')
     error(message('bioinfo:ngsbrowser:JavaMissing'));
end

if(nargin == 0) %start browser with no alignment
    openBrowser;
    return;
end

bioMapVarName = inputname(1);
    
% Validate BioMap object:
if ~(isobject(bioMapVar) && isa(bioMapVar, 'BioMap'))
    error(message('bioinfo:ngsbrowser:NotAValidBioMapObject'));
elseif numel(bioMapVar)~=1
    error(message('bioinfo:ngsbrowser:NotScalarInput'));
elseif numel(bioMapVar.SequenceDictionary)>1
    error(message('bioinfo:ngsbrowser:BioMapDataMultipleReferences',bioMapVarName));
elseif bioMapVar.NSeqs==0
    error(message('bioinfo:ngsbrowser:BioMapDataEmpty',bioMapVarName));
else
    st = getStart(bioMapVar);
    if ~issorted(st(st>0))
        error(message('bioinfo:ngsbrowser:BioMapDataNotSorted',bioMapVarName));
    end
end
        
openBrowser(bioMapVarName);

end

function openBrowser(varargin)
% OPENGBROWSER pass control to ngsbrowser java classes
try
    com.mathworks.toolbox.bioinfo.bioinfoservices.GenomeBrowserServices.openBrowser(varargin{1:end});
catch ME 
    error(message('bioinfo:ngsbrowser:GBJavaError'));
end
end

