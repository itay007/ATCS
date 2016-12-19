function [b,matrixInfo] = pam(N,varargin)
%PAM returns members of the PAM family of scoring matrices.
%
%   B = PAM(N) returns the PAMN matrix. Supported values of N are
%   10:10:500. Default ordering of the outputs is
%   A R N D C Q E G H I L K M F P S T W Y V B Z X *.
%
%   B = PAM(...,'EXTENDED',false) returns the scoring matrix for the 20
%   amino acids only and not for the extended symbols B,Z,X, and *.
%
%   B = PAM(...,'ORDER',ORDER) returns the PAMN matrix ordered by
%   the amino acid sequence ORDER. If ORDER does not contain the extended
%   characters B,Z,X, and *, then these characters are not returned.
%
%   [B,MATRIXINFO] = PAM(N) returns a structure of information about
%   the PAMN matrix with fields Name, Scale, Entropy, Expected, and
%   Order.
%
%   Examples:
%
%      PAM50 = pam(50)
%      PAM250 = pam(250,'order','CSTPAGNDEQHRKMILVFYW')
%
%   See also BLOSUM, DAYHOFF, GONNET, NWALIGN, PAM250, SWALIGN.

%   Reference: Dayhoff M.O., Schwartz R. and Orcutt B.C. (1978) Atlas of
%   protein sequence and structure. Vol. 5, Suppl. 3, Ed. M. O. Dayhoff.
%

%   Copyright 2003-2012 The MathWorks, Inc.


extended=true;
ordered=false;
possibleN = 10:10:500;

names='ARNDCQEGHILKMFPSTWYVBZX*';
intOrder = aa2int(names);

if ~ismember(N,possibleN)
    error(message('bioinfo:pam:InvalidBlosumN', N));
end

if nargin > 1
    if rem(nargin,2)== 0
        error(message('bioinfo:pam:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'extended','order'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname))); 
        if isempty(k)
            error(message('bioinfo:pam:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:pam:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % extended
                    extended = bioinfoprivate.opttf(pval);
                    if isempty(extended)
                        error(message('bioinfo:pam:InputOptionNotLogical', upper( char( okargs( k ) ) )));
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
                        error(message('bioinfo:pam:ShortOrder'));
                    end
                    if any(intOrder == 0)
                        error(message('bioinfo:pam:InvalidOrder'));
                    end
            end
        end
    end
end

[b,matrixInfo] = feval(sprintf('pam%d',N));

% permute matrix if necessary.
if ordered || ~extended
    if extended == false
        intOrder(intOrder>20) = [];
    end
    b = b(intOrder,intOrder);
    matrixInfo.Order = int2aa(intOrder);
end
