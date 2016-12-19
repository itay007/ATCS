function out = rnaconvert(s)
%RNACONVERT transforms RNA secondary structure between bracket and matrix
%notation
%
%   RNACONVERT(S) transforms the secondary structure S into matrix notation
%   (if S is in bracket notation) or into bracket notation (if S is in matrix
%   notation).
%
%   Example 1
%   S1 = '(((..((((.......)))).((.....)).))).';
%   M1 = rnaconvert(S1);
%
%   Example 2
%   M2 = zeros(12);
%   M2(1,12) = 1; M2(2,11) = 1; M2(3,10) = 1; M2(4,9) = 1;
%   S2 = rnaconvert(M2);
%
%   See also RNADEMO, RNAFOLD, RNAPLOT.

%   Copyright 2007-2008 The MathWorks, Inc.


%--------------------------------------------------------------------------

if isempty(s)
    error(message('bioinfo:rnaconvert:EmptyInput'))
else
    if ischar(s)
        if all(s == '(' | s =='.' | s== ')')
            ternary = (s == '(')-(s == ')');
            %if ((ternary == 1) | (ternary == 0) | (ternary == -1))
                if (sum(ternary))
                    error(message('bioinfo:rnaconvert:InvalidInputUnbalanced'))
                elseif (any(cumsum(ternary) < 0))
                    error(message('bioinfo:rnaconvert:InvalidInputUnordered'))
                else
                    out = bracket2matrix(s);
                end
            %end
        else % symbols other than '(.)' are in the structure
            error(message('bioinfo:rnaconvert:InvalidInputBracket'))
        end
       
    elseif isnumeric(s)
        if isvector(s)
            error(message('bioinfo:rnaconvert:InvalidInputVector'))
        elseif (size(s,1) ~= size(s,2))
            error(message('bioinfo:rnaconvert:InvalidInputSquare'))
        elseif (any(any(s ~= triu(s))))
            error(message('bioinfo:rnaconvert:InvalidInputTriu'))
        else
            out = matrix2bracket(s);
        end
    else
        error(message('bioinfo:rnaconvert:UnknownInputFormat'));
    end
end

%--------------------------------------------------------------------------
% Subfunctions
%--------------------------------------------------------------------------

function bracket = matrix2bracket(matrix)
% Convert matrix into bracket notation.

bracket = repmat('.', 1, size(matrix,1));
[i,j] = find(matrix);
bracket(i)='(';
bracket(j)=')';

%--------------------------------------------------------------------------

function matrix = bracket2matrix(bracket)
% Convert bracket notation into a (directed) matrix.

N = numel(bracket);  % number of bases
matrix = zeros(N,N); % matrix(i,j)=1 iff i<->j
topstack = 0;
openstack = [];

for i=1:N
    if (bracket(i)=='(')
        topstack = topstack + 1;
        openstack(topstack) = i; %#ok<AGROW>
    elseif (bracket(i)==')')
        %matrix(i,openstack(topstack)) = 1;
        matrix(openstack(topstack),i)=1; 
        topstack = topstack-1;
    end
end
