function cpgisland =  cpgisland(sequence, varargin)
%CPGISLAND locates CpG islands in a DNA sequence.
%
%   CPGISLAND(SEQ) finds CpG islands by tallying moving average of 100 DNA
%   bases with GC content greater than 50% and CpGobserved/CpGexpected
%   ratio greater than 60%. The default minimum island size is 200 bases.
%   NOTE: CPGISLAND does not count ambiguity bases, or gaps
%
%   CPGISLAND(...,'WINDOW',WINDOW) sets the window size for calculating GC
%   content and CpGobserved/CpGexpected averages for a sequence. The
%   default is 100 bases.
%
%   CPGISLAND(...,'MINISLAND',ISLANDSIZE) determines the minimum island
%   size to report. The default is 200 bases.
%
%   CPGISLAND(...,'CPGOE',CPGOE) specifies the cutoff for the minimum
%   CpGobserved/CPGexpected ratio. This is defined as
%    CPGobserved/CpGexpected =
%        (number of CpG's * Length)/( number of G's * number of C's)
%   The default is 0.6.
%
%   CPGISLAND(...,'GCMIN',GCMIN) specifies the minimum GC content for
%   each window. The default is 0.5.
%
%   CPGISLAND(...,'PLOT',true) will plot GC content, CPGoe content, CpG
%   islands greater than the minimum island size, and all potential CpG
%   islands for the specified criteria.
%
%   Example:
%
%       S = getgenbank('AC156455')
%       cpgisland(S.Sequence,'PLOT',true)
%
%   See also BASECOUNT, NTDENSITY, SEQSHOWORFS.

%   References:
%   Gardiner-Garden, M and Frommer, M. 1987. CpG islands in
%   vertebrate genomes. J.Mol.Biol. 196: 261-282.

%   Copyright 2003-2012 The MathWorks, Inc.

%check seq format

if isstruct(sequence)
    sequence = upper(bioinfoprivate.seqfromstruct(sequence));
elseif(~bioinfoprivate.isnt(sequence))
    error(message('bioinfo:CPGIsland:IncorrectSequenceType'));
elseif(isnumeric(sequence))
    sequence = upper(int2nt(sequence));
else
    sequence = upper(sequence);
end
%defaults
window = 100;
min_island_size = 200;
cgoe_min = 0.6;
gc_min = 0.5;
plotcpg = false;
%%% check arguments
if nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:CPGIsland:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'window','minisland', 'cpgoe', 'gcmin', 'plot'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        if(isstruct(pval))
            error(message('bioinfo:CPGIsland:StructParamError'));
        end
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:CPGIsland:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:CPGIsland:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % window
                    if(pval > 1 && pval < length(sequence))
                        window =  pval;
                    else
                        error(message('bioinfo:CPGIsland:IncorrectWindowSize'));
                    end
                case 2  % minisland
                    if(pval > 0 && pval < length(sequence))
                        min_island_size =  pval;
                    else
                        error(message('bioinfo:CPGIsland:IncorrectIslandSize'));
                    end
                case 3  % cpgoe
                    if(pval > 0 && pval < 1)
                        cgoe_min = pval;
                    else
                        error(message('bioinfo:CPGIsland:IncorrectCPGOESize'));
                    end
                case 4 % gcmin
                    if(pval > 0 && pval < 1)
                        gc_min = pval;
                    else
                        error(message('bioinfo:CPGIsland:IncorrectGCSize'));
                    end
                case 5 % plot
                    plotcpg = pval;
            end
        end
    end
end

b = ones(1,window,'single');

% Larsen et. al. further smooth the window calculations
% LarsenWindow = ones(1,10,'single')./10;

win_g = conv(single(sequence == 'G'),b,'same'); % G count in each sliding window
win_c = conv(single(sequence == 'C'),b,'same'); % C count in each sliding window
win_gc = win_c+win_g; % GC count in each sliding window

% Larsen et. al. further smooth the window calculations
% win_gc = conv(win_gc,LarsenWindow,'same');

% Expected CpG count times window size in each sliding window
win_cpg_exp = win_g .* win_c;

% Observed CpG count in each sliding window
% win_cpg_obs = conv(single([(sequence(1:end-1) == 'C' & sequence(2:end) == 'G') false]),b,'same');

% Observed/Expected ration in each sliding window
win_cpg_oeratio = (single(window) .* conv(single([(sequence(1:end-1) == 'C' & sequence(2:end) == 'G') false]),b,'same')) ./ win_cpg_exp;
win_cpg_oeratio(isnan(win_cpg_oeratio))=0;

% Larsen et. al. further smooth the window calculations
% win_cpg_oeratio = conv(win_cpg_oeratio,LarsenWindow,'same');

% CpG rich regions
cpg_rich = (win_gc>=(gc_min*window)) & (win_cpg_oeratio>=cgoe_min);

% Find limits of all islands
starts = find(diff([single(0) cpg_rich])>0)';
stops = find(diff([cpg_rich single(0)])<0)';

% Filter out islands by min island size
i = (stops-starts)>=(min_island_size-1);
starts = starts(i);
stops = stops(i);


if(plotcpg)
    figure;
    subplot(4,1,1);
    plot(win_gc./window);
    title('GC content');
    
    subplot(4,1,2);
    plot(win_cpg_oeratio);
    title('CPGoe content');
    
    subplot(4,1,3);
    cpg_islands = zeros(size(cpg_rich),'single');
    cpg_islands(starts) = 1;
    cpg_islands(stops+1) = -1;
    cpg_islands = cumsum(cpg_islands)>0;
    plot(cpg_islands);
    title(sprintf('CpG islands > %d bases',min_island_size ));

    subplot(4,1,4);
    plot(cpg_rich);
    title('All CpG islands');
end

if(isempty(starts))
    cpgisland.Starts = [];
    cpgisland.Stops = [];
else
    cpgisland.Starts = starts';
    cpgisland.Stops  = stops';
end


