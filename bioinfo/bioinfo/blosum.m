function [b,matrixInfo] = blosum(N,varargin)
%BLOSUM return members of the BLOSUM scoring matrix family.
%
%   B = BLOSUM(N) returns the BLOSUMN matrix. Supported values of N are
%   30:5:90,62 and 100. Default ordering of the outputs is
%   A R N D C Q E G H I L K M F P S T W Y V B Z X *.
%
%   B = BLOSUM(...,'EXTENDED',false) returns the scoring matrix for the 20
%   amino acids only and not for the extended symbols B,Z,X, and *.
%
%   B = BLOSUM(...,'ORDER',ORDER) returns the Blosum N matrix ordered by
%   the amino acid sequence ORDER. If ORDER does not contain the extended
%   characters B,Z,X, and *, then these characters are not returned.
%
%   [B,MATRIXINFO] = BLOSUM(N) returns a structure of information about
%   the BLOSUM N matrix with fields Name, Scale, Entropy, Expected, and
%   Order.
%
%   Examples:
%
%      B50 = blosum(50)
%      B75 = blosum(75,'order','CSTPAGNDEQHRKMILVFYW')
%
%   See also BLOSUM62, DAYHOFF, GONNET, NWALIGN, PAM, SWALIGN.

%   Reference: Henikoff S. and Henikoff J.G. (1992). Amino acid
%   substitution matrices from protein blocks. Proc. Natl. Acad. Sci. USA
%   89:10915-10519

%   Copyright 2002-2012 The MathWorks, Inc.


extended=true;
ordered=false;
possibleN=[30:5:90, 62, 100];

names='ARNDCQEGHILKMFPSTWYVBZX*';
intOrder = aa2int(names);

if ~ismember(N,possibleN)
    error(message('bioinfo:blosum:InvalidBlosumN', N));
end

if nargin > 1
    if rem(nargin,2)== 0
        error(message('bioinfo:blosum:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'extended','order'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname))); 
        if isempty(k)
            error(message('bioinfo:blosum:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:blosum:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % extended
                    extended = bioinfoprivate.opttf(pval);
                    if isempty(extended)
                        error(message('bioinfo:blosum:InputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 2 %order
                    order=pval;
                    ordered=true;
                    if ischar(order)
                        intOrder = aa2int(order);
                    else
                        intOrder = order;
                    end

                    if numel(intOrder)<20
                        error(message('bioinfo:blosum:ShortOrder'));
                    end
                    if any(intOrder == 0)
                        error(message('bioinfo:blosum:InvalidOrder'));
                    end
            end
        end
    end
end

[b,matrixInfo] = feval(sprintf('blosum%d',N));

% permute matrix if necessary.
if ordered || ~extended
    if extended == false
        intOrder(intOrder>20) = [];
    end
    b = b(intOrder,intOrder);
    matrixInfo.Order = int2aa(intOrder);
end
