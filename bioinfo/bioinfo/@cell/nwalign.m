function varargout = nwalign(seqs1,seqs2,varargin)
%NWALIGN performs Needleman-Wunsch global alignment of two sequences.
%
%   Implementation of NWALIGN for cell arrays. See NWALIGN for more info.

% Copyright 2005-2009 The MathWorks, Inc.


% check that none of the parameter-value input arguments are cells,
% otherwise we would iterate over the same function
if any(cellfun('isclass',varargin,'cell'))
    error(message('bioinfo:cell:nwalign:cellInputArguments'))
end

if nargout>2
    error(message('bioinfo:cell:nwalign:TooManyOutputs'))
end

if ~iscell(seqs1)
    seqs1 = {seqs1};
end
nseqs1 = numel(seqs1);

if isnumeric(seqs2) && isscalar(seqs2)
    % this form is used to compute N columns of a half-matrix full pairwise
    % alignment, output is a single column vector
    if ~any(1:nseqs1-1 == seqs2)
        error(message('bioinfo:cell:nwalign:InvalidNumberOfColumns'))
    end
    ncols = seqs2;
    k = 1;
    sc = zeros(ncols/2*(2*nseqs1-ncols-1),1);
    switch nargout
        case {0,1}
            for i = 1:ncols
                for j = i+1:nseqs1
                    sc(k) = nwalign(seqs1{i},seqs1{j},varargin{:});
                    k = k+1;
                end
            end
            varargout = {sc};
        case 2
            al = cell(ncols/2*(2*nseqs1-ncols-1),1);
            for i = 1:ncols
                for j = i+1:nseqs1
                    [sc(k), al{k}] = nwalign(seqs1{i},seqs1{j},varargin{:});
                    k = k+1;
                end
            end
            varargout = {sc al};
    end
else
    % this form is used to align all the sequences from seqs1 to all the
    % sequences from seqs2, the output is nseqs1 x nseqs2, and sequences should
    % not coexist in both sets
    if ~iscell(seqs2)
        seqs2 = {seqs2};
    end

    nseqs2 = numel(seqs2);

    sc = zeros(nseqs1,nseqs2);
    switch nargout
        case {0,1}
            for i = 1:nseqs1
                for j = 1:nseqs2
                    sc(i,j) = nwalign(seqs1{i},seqs2{j},varargin{:});
                end
            end
            varargout = {sc};
        case 2
            al = cell(nseqs1,nseqs2);
            for i = 1:nseqs1
                for j = 1:nseqs2
                    [sc(i,j), al{i,j}] = nwalign(seqs1{i},seqs2{j},varargin{:});
                end
            end
            varargout = {sc al};
    end
end
