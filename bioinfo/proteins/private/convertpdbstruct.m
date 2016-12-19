function pdbstruct = convertpdbstruct(oldpdbstruct, varargin)
%CONVERTPDBSTRUCT converts an old MATLAB PDB structure to a new PDB
% structure.
% 
%   PDBOUT = CONVERTPDBSTRUCT(OLDPDB) converts an old MATLAB PDB structure
%   OLDPDB to a new PDB structure PDBOUT. PDB structures generated from
%   PDBREAD and GETPDB have been changed since Bioinformatics Toolbox
%   Version 2.5. See PDBREAD and GETPDB doc for more information.
%
%   PDBOUT = CONVERTPDBSTRUCT(OLDPDB, FILENAME) warns about the structure
%   change used by the Matlab file FILENAME.
% 
%   See also GETPDB,PDBDISTPLOT, PDBREAD. 

%   Copyright 2007-2012 The MathWorks, Inc.


% Check oldpdbstruct is a structure.
if ~isstruct(oldpdbstruct)
    error(message('bioinfo:convertpdbstruct:InvalidInputStructure'))
end

filename = [];
if nargin == 2
    filename = varargin{1};
    
    % Check filename is string.
    if ~ischar(filename) || isempty(filename),
        error(message('bioinfo:convertpdbstruct:InvalidFileName'));
    end
end

% Field names for the coordinate information.
cofields = {'Atom',...
            'AtomSD',...
            'AnisotropicTemp',...
            'AnisotropicTempSD',...
            'Terminal',...
            'HeterogenAtom'};
        
realfields = isfield(oldpdbstruct, cofields); 
pdbstruct = oldpdbstruct;

if any(realfields)
    if isfield(pdbstruct, 'Model')
        N = numel(pdbstruct.Model);
        pdbstruct.Model = [];
        pdbstruct.Model.ModelNum = N;
    end
    
    for i=1:numel(cofields)
        if realfields(i)
            pdbstruct.Model.(cofields{i}) = oldpdbstruct.(cofields{i});
        end
    end
    
    pdbstruct = rmfield(pdbstruct, cofields(realfields));
    if ~isempty(filename)
        warning(message('bioinfo:convertpdbstruct:OldPDBStructureDetected', filename));
    end
end
