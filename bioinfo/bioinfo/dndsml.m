function [dn ds like params] = dndsml(s1,s2,varargin)
%DNDSML (non)-synonymous substitution rates by the maximum likelihood method.
%
%   [DN DS LIKE] = DNDSML(SEQ1,SEQ2) estimates synonymous and nonsynonymous
%   substitution rates between two homologous sequences by the maximum
%   likelihood method. This function returns the nonsynonymous substitution
%   rate (DN), the synonymous substitution rate (DS), and the likelihood
%   of this estimate (LIKE).
%
%   DNDSML implements the Goldman and Yang method '94. It estimates (using
%   maximum likelihood) an explicit model for codon substitution that takes
%   into account transition/transversion rate bias and base/codon
%   frequency bias. Then it uses the model to correct synonymous and
%   non-synonymous counts to account for multiple substitutions at the
%   same site. The maximum likelihood method is best suited when the sample
%   size is significant and when the sequences being compared may have
%   transition/transversion rate biases and base/codon frequency biases.
%   Any codons that include gaps or ambiguous nucleotide characters are
%   excluded from calculation. This analysis considers the number of codons
%   in the shortest sequence and assumes sequences do not have frame shifts.
%
%   DNDSML(...,'GENETICCODE',ID) calculates synonymous and nonsynonymous
%   substitution rates using the specified genetic code integer ID. The
%   default is Standard or 1.
%
%   DNDSML(...,'VERBOSE',true) turns verbosity on, use this option to
%   verify that the input sequences are codon-aligned. Default is false.
%
%   Examples:
%
%     % Estimate the substitution rates between the gag genes of two HIV
%     % strains:
%     gag1 = getgenbank('L11768')
%     gag2 = getgenbank('L11770')
%     [dn ds like] = dndsml(gag1, gag2)
%
%     % Extract the coding region of the neuraminidase protein (NA) of two
%     % strains of the Influenza A virus (H5N1), align them as aminoacids
%     % and copy the gaps to the nucleotide sequences for calculating the
%     % synonymous and nonsynonymous substitutions rates:
%     hk01 = getgenbank('AF509094');
%     vt04 = getgenbank('DQ094287');
%     hk01_cds = featuresparse(hk01,'feature','CDS','Sequence',true);
%     vt04_cds = featuresparse(vt04,'feature','CDS','Sequence',true);
%     [sc,al]=nwalign(nt2aa(hk01_cds),nt2aa(vt04_cds),'extendgap',1);
%     hk01_aligned = seqinsertgaps(hk01_cds,al(1,:))
%     vt04_aligned = seqinsertgaps(vt04_cds,al(3,:))
%     [dn,ds] = dndsml(hk01_aligned,vt04_aligned,'verbose',true)
%
%   See also DNDS, GENETICCODE, NT2AA, SEQINSERTGAPS, SEQPDIST.

%References:
%[1]  Yang, Z and Nielsen, R. "Estimating Synonymous and Nonsynonymous
%     Substitution Rates Under Realistic Evolutionary Models."
%     Mol.Biol.Evol. 17(1):32-43 (2000).
%[2]  Goldman, N and Yang, Z. "A Codon-based Model of Nucleotide
%     Substitution for Protein-coding DNA Sequences." Mol.Biol.Evol.
%     11(5):725-736 (1994).

%   Copyright 2003-2012 The MathWorks, Inc.



% default args
gene_code =1;
verbose = false;

if isstruct(s1)
    s1 =bioinfoprivate.seqfromstruct(s1);
end
if isstruct(s2)
    s2 = bioinfoprivate.seqfromstruct(s2);
end

if(~(bioinfoprivate.isnt(s1) && bioinfoprivate.isnt(s2)))
    error(message('bioinfo:dndsml:IncorrectSequenceType'));
end

if ~isnumeric(s1)
    s1 = nt2int(s1);
end
if ~isnumeric(s2)
    s2 = nt2int(s2);
end

% use common shorter length multiple of 3
num_codons =  min(floor(numel(s1)/3),floor(numel(s2)/3));
num_nt = num_codons * 3;

s1 = reshape(s1(1:num_nt),3,num_codons);
s2 = reshape(s2(1:num_nt),3,num_codons);

% gaps and ambiguous symbols are ignored, remove the whole codon
h = all(s1>=1,1) & all(s1<=4,1) & all(s2>=1,1) & all(s2<=4,1);
s1 = s1(:,h);
s2 = s2(:,h);

if num_codons<1
    error(message('bioinfo:dndsml:SequenceTooShort'));
end
%check varargin
if nargin > 2
    if rem(nargin,2) == 1
        error(message('bioinfo:dndsml:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'GENETICCODE','VERBOSE'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:dndsml:UnknownParameterName', pname));
        else
            switch(k)
                case 1  % GENETICCODE
                    if (isnumeric(pval) && isreal(pval))
                        gene_code =  pval;
                    else
                        error(message('bioinfo:dndsml:BadGeneCodeParam'));
                    end
                case 2 % VERBOSE
                    verbose = bioinfoprivate.opttf(pval);
                    if isempty(verbose)
                        error(message('bioinfo:dndsml:verboseInputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
            end
        end
    end
end

seq1 = int2nt(s1(:)');
seq2 = int2nt(s2(:)');
aaseq1 = nt2aa(seq1,'GeneticCode',gene_code);
aaseq2 = nt2aa(seq2,'GeneticCode',gene_code);

% remove stop codons from sequences
h = find(aaseq1=='*' | aaseq2=='*');
if ~isempty(h)
    aaseq1(h)=[];
    aaseq2(h)=[];
    seq1([(h-1)*3+1 (h-1)*3+2 (h-1)*3+3]')=[];
    seq2([(h-1)*3+1 (h-1)*3+2 (h-1)*3+3]')=[];
end
num_codons = numel(aaseq1);
num_nt = num_codons * 3;

if verbose
    disp('DNDSML: ')
    disp('Codons considered in the computations:')
    disp(seq1)
    disp(seq2)
    disp('Translations:')
    disp(sprintf('%c  ',aaseq1))
    disp(sprintf('%c  ',aaseq2))
end

%[pi codons intcodon] = get_eq_codon_freq(seq1, gene_code);
[pi codons intcodon] = get_eq_codon_freq([seq1 seq2], gene_code);
pi = repmat(pi',length(pi),1);
aa = nt2aa(cell2mat(codons'),'ALTERNATIVESTARTCODONS',false,'GENETICCODE',gene_code); % translation
sites = get_sites(seq1,seq2,codons,num_nt); % get sites for likelihood calculation
indQ  = build_indQ(intcodon,aa);  %   stv: 1, sts: 2, nstv: 3, nsts: 4

% get initial estimates
[dn ds] = dnds(seq1,seq2,'adjuststops',true); % get initial guess of omega

if isnan(dn) || isnan(ds)
    like = NaN;
    params = [NaN NaN NaN];
    return
end

[kappa, t] = HKY85_F84(seq1,seq2,num_nt); % get initial guess at kappa and t

if verbose
    disp(sprintf('\nInitial estimates: Kappa=%f, dn=%f, ds=%f, t=%f',kappa,dn,ds,t))
end

if isnan(kappa) % a warning was already thrown inside HKY85_F84, this case happens when sequences are too divergent
    ds = NaN;
    dn = NaN;
    like = NaN;
    params = [NaN NaN t];
    return
end

% Finds the max likelihood
if dn==0 || ds==0 % handle edge case with absence of synonynous counts or nonsynomous counts
    if isinf(kappa) || kappa==0 % handle case with also edge tranversion/transition rate ratios
        if dn==0 && ds==0 % equal sequences, just figure out the likelihood value
             params = [0 0 0];
             like = get_likelihood(1,1,0,sites,pi,indQ);
             exitFlag = 1;
        else % estimate only t by ML
            [params like exitFlag] = fminsearch(@(x)get_likelihood(kappa,abs(dn/ds),x,sites,pi,indQ),t,optimset('Display','off'));
            params = [kappa abs(dn/ds) params];
        end
    else  % estimate only kappa and t by ML
        [params like exitFlag] = fminsearch(@(x)get_likelihood(x(1),abs(dn/ds),x(2),sites,pi,indQ),[kappa t],optimset('Display','off'));
        params = [params(1) abs(dn/ds) params(2)];
    end
elseif isinf(kappa) || kappa==0 % handle case with edge tranversion/transition rate ratios
    [params like exitFlag] = fminsearch(@(x)get_likelihood(kappa,x(1),x(2),sites,pi,indQ),[dn/ds t],optimset('Display','off'));
    params = [kappa params(1) params(2)];
else % estimate kappa, omega and t by ML
   [params like exitFlag] = fminsearch(@(x)get_likelihood(x(1),x(2),x(3),sites,pi,indQ),[kappa dn/ds t],optimset('Display','off'));
end

if exitFlag == 0
   warning(message('bioinfo:dndsml:noConverging'))
   [params like] = fminsearch(@(x)get_likelihood(abs(kappa),abs(dn/ds),x,sites,pi,indQ),t);
   params = [kappa abs(dn/ds) params];
end

if verbose
   disp(sprintf('ML estimates: Kappa=%f, omega(dn/ds)=%f, t=%f',params))
end

% Calculates dn ds
Q = build_Q(params(1),params(2),pi,indQ);
[rhos_star rhon_star] = get_rho(Q,aa,pi);
Q = build_Q(params(1),1,pi,indQ);
[rhos_inf rhon_inf] = get_rho(Q,aa,pi);
if rhos_star==0
    ds=0;
else
    ds = (params(3) * rhos_star)/(3 * rhos_inf);
end
if rhon_star==0
    dn = 0;
else
    dn = (params(3) * rhon_star)/(3 * rhon_inf);
end
like = -like;

function Q = build_Q(kappa,omega,pi,Q)
% Returns the codon substitution rate matrix
if isinf(omega)
    K = 0; omega = 1;
else
    K = 1;
end
if isinf(kappa)
    J = 0; kappa = 1;
else
    J = 1;
end
syntv = (Q == 1);
synts = (Q == 2);
nsyntv = (Q == 3);
nsynts = (Q == 4);
Q(syntv) = J*K*pi(syntv);                % (1) Yang&Nielsen, 2000
Q(synts) = K*kappa*pi(synts);            % (1) Yang&Nielsen, 2000
Q(nsyntv) = J*omega*pi(nsyntv);          % (1) Yang&Nielsen, 2000
Q(nsynts) = omega*kappa*pi(nsynts);      % (1) Yang&Nielsen, 2000
Q = Q + diag(-sum(Q,2));                 % (2) Yang&Nielsen, 2000
D = -sum(diag(pi).*diag(Q));
if D~=0
    Q = Q ./ D ;        % (3) Yang&Nielsen, 2000
else
    Q(:) = 0;
end
    
% sum(Q,2) == 0 
% - sum(diag(pi).*diag(Q)) == diag(pi)'*sum(tril(Q,-1)+triu(Q,1),2)  == 1

function like = get_likelihood(kappa,omega,t,n,pi,Q)
Q = build_Q(kappa,omega,pi,Q);
P = expm(Q*t);                           % (4) Yang&Nielsen, 2000
F = pi.*P;                               % (5) Yang&Nielsen, 2000
Fn = max(0,F(n>0)); % removes small negative numbers that may appear due to 
                    % rounding errors, by setting these to 0, log(0)->-Inf 
                    % and like->Inf avoiding that FMINSEARCH picks this as
                    % a viable solution.  
like = -sum((nonzeros(n).*log(Fn)));    % (6) Yang&Nielsen, 2000


function indQ  = build_indQ(intcodon,aa)
% Returns codon pair types corresponding to the matrix Q as follows:
%   - synonymous traversions: 1
%   - synonymous transitions: 2
%   - nonsynonymous traversions: 3
%   - nonsynonymous transitions: 4
indQ = zeros(length(intcodon),length(intcodon));
for i = 1:length(intcodon)
    diffs = [intcodon(i,1)~=intcodon(:,1) ...
        intcodon(i,2)~=intcodon(:,2) ...
        intcodon(i,3)~=intcodon(:,3)];
    trans = [rem(intcodon(i,1),2)~=rem(intcodon(:,1),2) ...
        rem(intcodon(i,2),2)~=rem(intcodon(:,2),2) ...
        rem(intcodon(i,3),2)~=rem(intcodon(:,3),2)];
    ind = sum(diffs,2) == 1;
    syn = aa(i) == aa;
    indQ(i,any(trans&diffs,2)&ind&syn') = 1;
    indQ(i,~any(trans&diffs,2)&ind&syn') = 2;
    indQ(i,any(trans&diffs,2)&ind&~syn') = 3;
    indQ(i,~any(trans&diffs,2)&ind&~syn') = 4;
end

function [rhos rhon] = get_rho(Q,aa,pi)
% Returns the proportions of synonymous and nonsynonymous substitutions
% (7) Yang&Nielsen, 2000
nsaa = repmat(aa,numel(aa),1)~=repmat(aa',1,numel(aa));
rhon = sum(diag(pi)'*(Q.*nsaa));
% fixing degenerative rounding errors in edge cases
if rhon<(61*eps)
    rhon = 0; rhos = 1;
elseif (1-rhon)<(61*eps)
    rhon = 1; rhos = 0;
else
    rhos = 1-rhon;
end
% saa = bsxfun(@(a1,a2) a1==a2, aa, aa');
% rhos == sum(diag(pi)'*((tril(Q,-1)+triu(Q,1)).*saa))

function sites = get_sites(seq1, seq2, codons, num_nt)
% Returns a matrix with site counts for each codon pair ignoring gaps
sites = zeros(length(codons), length(codons));
for i = 1:3:num_nt-2
    ind1 = strcmp(seq1(i:i+2),codons);
    ind2 = strcmp(seq2(i:i+2),codons);
    if(any(ind1) && any(ind2))
        sites(ind1,ind2) = sites(ind1,ind2) +1;
    end
end

function [f codons int_codon] = get_eq_codon_freq(seq,gene_code)
% Calculates equilibrium codon frequency via product of nucleotide
% frequency and adjust frequencies for terminator codons
[lastwmsg,lastwid]=lastwarn;
warnState=warning('off'); %#ok
[codon_freq cod_array] = codoncount(seq);
warning(warnState);
lastwarn(lastwmsg,lastwid);
nt_freq = ...
    [sum(sum(cod_array(1,:,:)))  sum(sum(cod_array(:,1,:))) sum(sum(cod_array(:,:,1))); ...
    sum(sum(cod_array(2,:,:)))  sum(sum(cod_array(:,2,:))) sum(sum(cod_array(:,:,2))); ...
    sum(sum(cod_array(3,:,:)))  sum(sum(cod_array(:,3,:))) sum(sum(cod_array(:,:,3))); ...
    sum(sum(cod_array(4,:,:)))  sum(sum(cod_array(:,4,:))) sum(sum(cod_array(:,:,4)))];
nt_freq  = nt_freq / sum(cod_array(:));
% get probs for termination codons
codon_table = geneticcode(gene_code);
codons = fieldnames(codon_table);
not_codon = (cellfun('length',codons) ~= 3);
codons(not_codon) = [];
translation = struct2cell(codon_table);
translation(not_codon) = [];
stpind = find(strcmp('*',translation));
stps = codons(stpind);
int_codon  = nt2int(cell2mat(codons));
int_stps  = nt2int(cell2mat(stps));
stp_total = sum(prod([nt_freq(int_stps(:,1),1) nt_freq(int_stps(:,2),2) nt_freq(int_stps(:,3),3)],2));
f = prod([nt_freq(int_codon(:,1),1) nt_freq(int_codon(:,2),2) nt_freq(int_codon(:,3),3)],2) ./ (1 - stp_total);
f(stpind) = [];
codons(stpind) = [];
int_codon(stpind,:) = [];

function [K_HKY85, t] = HKY85_F84(s1,s2,num_nt)
% Hasegawa distance taken from seqpdist modified according to Yang Nielsen
%
% Note: nucleotide frequencies are averaged between s1 and s2
%
% Let P and Q be the transition and transversional difference proportions
% and g the nucleotide frequencies over all the sequences:
%
% gr = g(1)+g(3);  % purines frequency
% gy = g(2)+g(4);  % pyrimidines frequency
% pr = g(1)*g(3);  % purines product
% py = g(2)*g(4);  % pyrimidines product
%
% Reference: Tamura and  Nei (1993) Molecular Biology and Evolution
% Reference: Yang and Nielsen(2000) Molecular Biology and Evolution

h = (s1(1:num_nt) ~='-') & (s2(1:num_nt)~='-');
map('ACGT- ') = 1:6;
temp1 = upper(char(s1(1:num_nt)));
temp2 = upper(char(s2(1:num_nt)));
g1 = accumarray(map(temp1(:))',1,[4 1]);
g2 = accumarray(map(temp2(:))',1,[4 1]);
g = (g1(1:4)+g2(1:4))/2;
g = g/ sum(g);
X  = accumarray(double([nt2int(s1(h))',nt2int(s2(h))']),1,[4 4]);
numPairs = sum(X(:));
numAGTransitions = sum(X([3 9]));
numCTTransitions = sum(X([8 14]));
numTransitions = numCTTransitions + numAGTransitions;
numTransversions = numPairs - sum(diag(X)) - numTransitions;
P = numTransitions / numPairs;
Q = numTransversions / numPairs;
gr = g(1)+g(3); % purines frequency
gy = g(2)+g(4); % pyrimidines frequency
pr = g(1)*g(3); % purines product
py = g(2)*g(4); % pyrimidines product

% Yang&Nielsen (8)
A = ( py+pr + (py*gr/gy+pr*gy/gr) * (1-(Q/(2*gy*gr))) - P/2 ) / (py/gy+pr/gr);
B = 1-Q/(2*gy*gr);

if B<=0
    warning(message('bioinfo:dndsml:QLarge'))
    K_HKY85 = NaN;
    t = NaN;
    return
elseif A<=0
    warning(message('bioinfo:dndsml:PLarge'))
    K_HKY85 = NaN;
    t = NaN;
    return
end

% Yang&Nielsen (9) (A and B must be positive !)
a = -log(A);
b = -log(B);

if Q==0 %(b==0)
    K_HKY85 = inf;
    t = 4*a*(py/gy+pr/gr); % taking the limit of Yang&Nielsen (10)
else
    K_F84 = a/b-1; % Yang&Nielsen (10)
    K_HKY85 = 1 + (((py/gy + pr/gr)*K_F84) / (py + pr));       % Yang&Nielsen (11)
    t = (4*py*(1+K_F84/gy) + 4*pr*(1+K_F84/gr) + 4*gy*gr) * b; % Yang&Nielsen (10)
end
K_HKY85 = max(0,K_HKY85);  % Kappa can't be negative, fixing rounding errors
t = max(0,t);              % t can't be negative, fixing rounding errors

