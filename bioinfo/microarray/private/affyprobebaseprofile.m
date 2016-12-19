function [A, stats, x] = affyprobebaseprofile(seqMatrix, intensities, varargin)
%AFFYPROBEBASEPROFILE estimates Affymetrix GeneChip probe base profiles
% from their sequences and probe intensities.
% 
%   A = AFFYPROBEBASEPROFILE(SEQMATRIX, INTENSITIES) estimates affinity
%   coefficients A using multiple linear regression. For a given type of
%   base, the positional effect is modeled as a polynomial of degree 3.
%   SEQMATRIX is an Nx25 matrix of the probe sequence bases in number, with
%   rows for the probes in a chip and every column is the position of 25.
%   INTENSITIES is the observed probe intensities. 
% 
%   [A, STATS] = AFFYPROBEBASEPROFILE(...) returns a vector of regression
%   statistics STATS containing, in the following order, the R-square
%   statistic, the F statistic and p-value, and an estimate of the error
%   variance.
%   
%   AFFYPROBEBASEPROFILE(...,'SHOWPLOT',TF) displays the affinities for
%   each base vs the positions if TF is true. Default is false.
% 
%   Example:
%       load affy95a_spikein
%       A = affyprobebaseprofile(smatrix, mmmatrix(:,1), 'showplot', true);
% 
%   See also AFFYPROBEAFFINITIES, AFFYPROBESEQREAD, CELINTENSITYREAD,
%   GCRMA, GCRMABACKADJ, REGRESS.

%   Copyright 2006 The MathWorks, Inc.


%   References:
%       [1] Naef F, Magnasco MO. "Solving the Riddle of the Bright
%       Mismatches: Labeling and Effective Binding in Oligonucleotide
%       Arrays" Physical Review E 68, 011906 (2003). 

if nargin < 2
    error(message('bioinfo:affyprobebaseprofile:NotEnoughInputs'));
end

showplot = false;

if ~isnumeric(intensities) || ~isreal(intensities) || ~isvector(intensities)
   error(message('bioinfo:affyprobebaseprofile:ProbeIntensityNotNumericAndRealVector')) 
end

if ~isa(seqMatrix, 'uint8')
    error(message('bioinfo:affyprobebaseprofile:SequenceMatrixNotUINT8')) 
end

N = numel(intensities);

if size(seqMatrix,1) ~= N || size(seqMatrix,2) ~= 25
    error(message('bioinfo:affyprobebaseprofile:SequenceMatrixWrongSize', N)); 
end

% Columnize the input if is a vector
intensities = intensities(:);

% Deal with the various input
if nargin > 2
    if rem(nargin, 2) == 1
        error(message('bioinfo:affyprobebaseprofile:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'showplot'};
    for j = 1:2:nargin-3
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs, numel(pname)));

        if isempty(k)
            error(message('bioinfo:affyprobebaseprofile:UnknownParameterName', pname));
        elseif length(k) > 1
            error(message('bioinfo:affyprobebaseprofile:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % showplot
                    showplot = bioinfoprivate.opttf(pval);
                    if isempty(showplot)
                        error(message('bioinfo:affyprobebaseprofile:InputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
            end % switch
        end % if
    end % for loop
end % if

% %------------------fitting to get affinity coefficients-----------------
Y = log2(intensities);
x = getprobebasematrix(seqMatrix);

% Get X matrix for regression of a polynomials of order 3
X = [ones(N,1) x];
[A, BINT, R, RINT, stats] = regress(Y, X);

% Plot
if showplot
    plotbaseprofiles(A(2:end))
end 

%-------------------Helper functions-------------------------------------%
function X = getprobebasematrix(seqMatrix)
%getprobebasematrix returns a basis matrix X of Nx12 from sequence profile
% Matrix. revX is the same matrix only the 13 is the reverse complement
%   N - number of probes. For a given type of base, its positional effect
%   is modeled as a polynomial of degree 3. We are fitting A, C, G bases
%   with experiment data, there are 12 parameters. T is zero.

% Get 3rd order polynormial base matrix - 25x4
k = (1:25)';
b = [ones(size(k)) k k.^2 k.^3];

X = zeros(size(seqMatrix,1), 12);
for i = 1:3
    seqX = (seqMatrix==i);
    X(:, 1+(i-1)*4:i*4) = seqX * b;
end


%------------------------------------------------------------------
function plotbaseprofiles(A)
% A - polynomial to 3-order coefficients for A, C, G
res = zeros(4, 25);
positions = 1:25;
colors = {'r', 'b', 'g', 'm'};

for i = 1:3
    a = A(1+(i-1)*4:i*4);
    a = flipud(a);
    rt = polyval(a, positions);
    res(i,:) = rt;
end

sumcol=sum(res,1)/4;
sumcol = repmat(sumcol, 4, 1);
res=res-sumcol;

affis = [min(min(res,[],2))-0.05, max(max(res, [],2))+0.05];

for base = 1:4
    hold on
    text(positions, res(base,:), int2nt(base), 'Color', colors{base})
end
haxis = gca;
set(haxis, 'xlim', [0 26], 'ylim', affis)
title('Position-dependent Affinity Base Profile')
xlabel('Position')
ylabel('Affinity')
hold off

% % %------------------------------------------------%
% % function A=initfitbaseprofile(Y, x) %#ok
% % % Naive fit of the base profile with 4 (bases) x 25 (positions) = 100 variables
% % 
% % A = zeros(1,100);
% % for i = 1:4
% %     X = [ones(numel(Y),1), x(:, 1+25*(i-1):25*i)];
% %     a = X\Y;
% %     A(1+25*(i-1):25*i) = a(2:26);
% % end
% % 
% % % plot the probe base profile
% % figure
% % colors={'r', 'b', 'g', 'm'};
% % positions = 1:25;
% % % for fitting 100
% % for base = 1:4
% %     hold on
% %     start = 1 + (base-1) * 25;
% %     stop = start + 24;
% % %     plot(positions, A(start:stop), [colors{base} 'x'])
% %     text(positions, A(start:stop), int2nt(base), 'Color', colors{base})
% % end
% % haxis = gca;
% % set(haxis, 'xlim', [0 26], 'ylim', [min(A)-0.1, max(A)+0.1])
% % hold off
% % 
% % 
