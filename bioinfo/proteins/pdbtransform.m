function pdbtx = pdbtransform(PDB, transf, varargin)
%PDBTRANSFORM applies a linear transformation to the 3D structure of a molecule.
%
%   PDBTRANSFORM(PDB, TRANSF) applies the linear transformation specified
%   in TRANSF to the coordinates of the molecule represented in PDB. PDB
%   can be a valid identifier in the PDB database, the name of a variable
%   containing a PDB-formatted MATLAB structure, or the name of a
%   PDB-formatted file. TRANSF is a MATLAB structure representing a linear
%   transformation, such as that returned by the functions PROCRUSTES and
%   PDBSUPERPOSE. TRANSF must have the following fields: 
%       T : orthogonal, rotation and reflection component
%       b : scale component 
%       c : translation component
%   
%   PDBTX = PDBTRANSFORM(PDB, TRANSF) returns the transformed PDB-formatted
%   structure in PDBTX.
% 
%   PDBTRANSFORM(..., 'SEGMENT', SEGMENTVALUE) specifies the extent to
%   which the linear transformation is to be applied. Valid options for
%   SEGMENTVALUE are: 'all' (the transformation is applied to the entire
%   PDB input) or a string specifying the boundaries and the chain of the
%   segment to consider, such as 'start-stop:chain', (the transformation is
%   applied to the specified segment only). The boundaries can be omitted
%   to indicate the entire chain.  Default is 'all'.
%
%   PDBTRANSFORM(..., 'MODEL', MODELVALUE) specifies the model to consider.
%   By default the first model is considered.
%
%   Example:
%
%   % Define a linear transformation 
%   transf.T = eye(3);  transf.b = 1;  transf.c = [11.8 -2.8 -32.3];
%   % Apply the linear transformation to chain B in thioredoxin structure (PDB ID 2trx)
%   pdbtx = pdbtransform('2trx', transf, 'segment', 'B');
%
%   See also GETPDB, MOLVIEWER, PDBREAD, PDBSUPERPOSE, PROCRUSTES.

%   Copyright 2008-2010 The MathWorks, Inc.


%=== input checking
bioinfochecknargin(nargin,2,mfilename);
[start, stop, chain, model, apply] = parse_inputs(varargin{:});

%=== check validity of the linear transformation specified
if isstruct(transf) && isfield(transf, 'T') && isfield(transf, 'b') && isfield(transf, 'c')
    if ~isequal(size(transf.T),[3 3]) 
        error(message('bioinfo:pdbtransform:InvalidComponentT'));
    elseif ~isscalar(transf.b)
        error(message('bioinfo:pdbtransform:InvalidComponentb'));
    elseif size(transf.c,2) ~= 3
        error(message('bioinfo:pdbtransform:InvalidComponentc'));
    end
else
    error(message('bioinfo:pdbtransform:InvalidTransformation'));
end

%=== read the pdb entry (structures, PDB files, or valid PDB IDs)
try 
    if (~isstruct(PDB))
        if exist(PDB,'file')
            pdbStruct = pdbread(PDB);
        else
            pdbStruct = getpdb(PDB);
        end
    else
        pdbStruct = convertpdbstruct(PDB);            
    end
catch theErr
    if ~isempty(strfind(theErr.identifier,'getpdb')) 
        rethrow(theErr);
    end
    error(message('bioinfo:pdbtransform:InvalidInput'));
end


%=== set default chain (first chain)
chainList = {pdbStruct.Sequence.ChainID}; % all chains

if isempty(chain)
    chain = chainList(1);
end

%=== check validity of selected chain
c = strmatch(chain, chainList, 'exact');

if isempty(c)
    error(message('bioinfo:pdbtransform:InvalidChain'));
end

%=== get specified model
modelStruct = pdbStruct.Model(model);

%=== initialize coords of transformed pdb
atomCoordsTx = [[modelStruct.Atom.X];  [modelStruct.Atom.Y];  [modelStruct.Atom.Z]]'; % coordinates of ALL atoms in the model


if apply == 1 % entire structure
    allAtomCoords = [[modelStruct.Atom.X]; [modelStruct.Atom.Y]; [modelStruct.Atom.Z]]';
    n = size(allAtomCoords, 1);
    atomCoordsTx = transf.b * allAtomCoords * transf.T + repmat(transf.c(1,:), n, 1);

else
    %=== determine range of residue positions in specified chain
    isChainAtom = [modelStruct.Atom.chainID]' == chain; % atoms in chain
    chainAtoms = modelStruct.Atom(isChainAtom); % data for atoms in chain
    
    if apply == 2 % specified chain

        chainAtomCoords = [[chainAtoms.X]; [chainAtoms.Y]; [chainAtoms.Z]]'; % nx3
        n = size(chainAtomCoords, 1);
        atomCoordsTx(isChainAtom,:) =  transf.b * chainAtomCoords * transf.T + repmat(transf.c(1,:), n, 1);
   
    else % apply  == 3, only segment
        %=== determine the residues in specified segment
                     
        res = [chainAtoms.resSeq]; % res positions for each atom in chain
        
        if start < min(res) || stop > max(res)
            error(message('bioinfo:pdbtransform:InvalidSegment'));
        end
        
        selectedRes = start:stop;
        isSelectedAtom = ismember(res, selectedRes); % atoms in selected res
        isSegAtom = isChainAtom;
        isSegAtom(isSegAtom) = isSelectedAtom; % logical w.r.t. the entire structure
        segAtoms = chainAtoms(isSelectedAtom);

        segAtomCoords = [[segAtoms.X]; [segAtoms.Y]; [segAtoms.Z]]'; % nx3
        n = size(segAtomCoords, 1);
        atomCoordsTx(isSegAtom,:) =  transf.b * segAtomCoords * transf.T + repmat(transf.c(1,:), n, 1);
    end
end


%=== write newly transformed coordinates
for i = 1:size(atomCoordsTx,1)
    modelStruct.Atom(i).X = atomCoordsTx(i,1);
    modelStruct.Atom(i).Y = atomCoordsTx(i,2);
    modelStruct.Atom(i).Z = atomCoordsTx(i,3);
end

pdbtx = pdbStruct;
pdbtx.Model(model) = modelStruct;

         
 %-------------------------------------------------------------------------
function [start, stop, chain, model, apply] = parse_inputs(varargin)
% Parse input PV pairs.
     
%=== Check for the right number of inputs
if rem(nargin,2)== 1
    error(message('bioinfo:pdbtransform:IncorrectNumberOfArguments', mfilename));
end

%=== Allowed inputs
okargs = {'model', 'segment'};

%=== Defaults
model = 1;   % model to consider
start = 0;   % start boundary of segment
stop = 0;    % stop boundary of segment
chain = '';  % chain to which segment belong
apply = 1;   % apply to the entire struct (1), the chain (2) or the segment (3)


for j=1:2:nargin
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1 % model
            if isnumeric(pval) &&  numel(pval) == 1
                model = pval;
            else
                error(message('bioinfo:pdbtransform:InvalidModelOption'));
            end
            
        case 2 % segment
                          
            if iscell(pval) 
                pval = pval{:}; 
            end
            
            if strcmpi(pval, 'all')
                apply = 1;
            else

                %=== validate syntax segment
                s = regexpi(pval, '^(\d+)-(\d+):(\w+)$', 'tokens');

                if ~isempty(s)
                    start = str2double(s{1}{1});
                    stop = str2double(s{1}{2});
                    chain = s{1}{3};
                    apply = 3;

                else % must be chain only
                    s = regexpi(pval, '^\w+$', 'match'); % assume chain are one alphanumeric char

                    if ~isempty(s)
                        chain = s{1};
                        apply = 2;
                    else
                        error(message('bioinfo:pdbtransform:InvalidSegmentOption'));
                    end

                end
                
                %=== validate boundaries
                if start > stop
                    error(message('bioinfo:pdbtransform:InvalidSegmentBoundaries'));
                end
            end
    end
end

