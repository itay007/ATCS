function [d_hat_n d_hat_s var_d_hat_n var_d_hat_s] = dnds(s1,s2,varargin)
%DNDS estimates synonymous and nonsynonymous substitutions rates.
%
%   [DN DS VARDN VARDS] = DNDS(SEQ1,SEQ2) estimates the synonymous and
%   nonsynonymous substitutions rates per site between two homologous
%   sequences by comparing codons using the Nei-Gojobori method, the
%   Li-Wu-Luo method or the Pamilo-Bianchi-Li method. This
%   function returns the nonsynonymous substitution rate (DN), the
%   synonymous substitution rate (DS), the variance for the nonsynonymous
%   substitutions rate (VARDN), and the variance for the synonymous
%   substitutions (VARDS). Any codon that includes gaps or ambiguous
%   nucleotide characters is excluded from calculation. This analysis
%   considers the number of codons in the shortest sequence and assumes
%   sequences do not have frame shifts.
%
%   DNDS(...,'GENETICCODE',ID) calculates synonymous and nonsynonymous
%   substitution rates using the specified genetic code integer ID. The
%   default is Standard or 1.
%
%   DNDS(...,'METHOD',method) allows synonymous/nonsynonymous substitution
%   rates to be calculated using one of the following approaches:
%     'NG'    - Nei-Gojobori method '86, uses on the number of synonymous
%   (default)   and non-synonymous substitutions and the number of
%               potentially synonymous and non-synonymous sites. Based on
%               the Jukes-Cantor model.
%     'LWL'   - Li-Wu-Luo method '85, uses the number of transitional and
%               transversional substitutions at three different levels of
%               degeneracy of the genetic code. Based on the Kimura's
%               2-parameter model.
%     'PBL'   - Pamilo-Bianchi-Li method '93, similar to LWL but with bias
%               correction. It should be used when the number of
%               transitions is much larger than the number of transversions.
%
%   DNDS(..., 'WINDOW',windowSize) performs calculations over a sliding
%   window of size specified by windowSize (in codons).
%
%   DNDS(...,'ADJUSTSTOPS',FALSE) treats stop codons as other codons in all
%   calculations. Default is TRUE, meaning that stop codons are excluded
%   from frequency tables and that paths containing stop codons are not
%   counted in the Nei-Gojobori method.
%
%   DNDS(...,'VERBOSE',true) turns verbosity on. Use this option to verify
%   that the input sequences are codon-aligned. Default is false.
%
%   Examples:
%
%    % Estimate the substitution rates between the gag genes of two HIV
%    % strains:
%    gag1 = getgenbank('L11768')
%    gag2 = getgenbank('L11770')
%    [dn ds vardn vards] = dnds(gag1,gag2)
%
%    % Extract the coding region of the neuraminidase protein (NA) of two
%    % strains of the Influenza A virus (H5N1), align them as aminoacids
%    % and copy the gaps to the nucleotide sequences for calculating the
%    % synonymous and nonsynonymous substitutions rates:
%    hk01 = getgenbank('AF509094');
%    vt04 = getgenbank('DQ094287');
%    hk01_cds = featuresparse(hk01,'feature','CDS','Sequence',true);
%    vt04_cds = featuresparse(vt04,'feature','CDS','Sequence',true);
%    [sc,al]=nwalign(nt2aa(hk01_cds),nt2aa(vt04_cds),'extendgap',1);
%    hk01_aligned = seqinsertgaps(hk01_cds,al(1,:))
%    vt04_aligned = seqinsertgaps(vt04_cds,al(3,:))
%    [dn,ds] = dnds(hk01_aligned,vt04_aligned,'verbose',true)
%
%   See also BIRDFLUDEMO, DNDSDEMO, DNDSML, GENETICCODE, NT2AA,
%   SEQINSERTGAPS, SEQPDIST.

%References:
%[1]  Nei, M and Gojobori, T. 1986. Simple Methods For Estimating The
%     Numbers Of Synonymous And Nonsynonymous Nucleotide Substitutions.
%     Mol.Biol.Evol. 3(5):418-426.
%[2]  Pamilo P. and Bianchi, N. 1993. Evolution Of The Zfx And Zfy Genes:
%     Rates and Interdependence Between The Genes. Mol.Biol.Evol.
%     10(2):271-281.
%[3]  Li, W, Wu, C, and Luo, C. 1985. A New Method For Estimating
%     Synonymous And Nonsynonymous Rates Of Nucleotide Substitution
%     Considering The Relative Likelihood Of Nucleotide And Codon Changes.
%     Mol.Biol.Evol.2(2):150-174.
%[4]  Nei, M and Jin, L. 1989. Variances Of The Average Numbers Of
%     Nucleotide Substitutions Within And Between Populations.
%     Mol.Biol.Evol. 6(3):290-300.
%[5]  Nei, M and Kumar, S. 2000. Molecular Evolution and Phylogenetics.
%     Oxford University Press. Chapter 4.

%   Copyright 2003-2012 The MathWorks, Inc.


%default args
gene_code = 1;
method = 'NG';
verbose = false;
windowSize = inf;
slidingWindow = true;
adjustStops = true;

if isstruct(s1)
    s1 =bioinfoprivate.seqfromstruct(s1);
end
if isstruct(s2)
    s2 = bioinfoprivate.seqfromstruct(s2);
end

if(~(bioinfoprivate.isnt(s1) && bioinfoprivate.isnt(s2)))
    error(message('bioinfo:dnds:IncorrectSequenceType'));
end

if ~isnumeric(s1)
    s1 = nt2int(s1);
end
if ~isnumeric(s2)
    s2 = nt2int(s2);
end

% use common shorter length multiple of 3
num_codons_original =  min(floor(numel(s1)/3),floor(numel(s2)/3));
num_nt = num_codons_original * 3;

s1 = reshape(s1(1:num_nt),3,num_codons_original);
s2 = reshape(s2(1:num_nt),3,num_codons_original);


%check varargin
if nargin > 2
    if rem(nargin,2) == 1
        error(message('bioinfo:dnds:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'GENETICCODE','METHOD','VERBOSE','WINDOW','ADJUSTSTOPS'};
    okmeth = {'NG', 'LWL', 'PBL'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:dnds:UnknownParameterName', pname));
        else
            switch(k)
                case 1  % GENETICCODE
                    try
                        geneticcode(pval);
                    catch allExceptions
                        error(message('bioinfo:dnds:BadGeneCodeParam'));
                    end
                    gene_code = pval;
                case 2  % METHOD
                    if (ischar(pval) && any(strcmpi(pval,okmeth)))
                        method = upper(pval);
                    else
                        error(message('bioinfo:dnds:BadMethodParam'));
                    end
                case 3 % VERBOSE
                    verbose = bioinfoprivate.opttf(pval);
                    if isempty(verbose)
                        error(message('bioinfo:dnds:verboseInputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 4 % WINDOW
                    windowSize = pval;
                    if ~isscalar(windowSize) || windowSize<=0 || islogical(windowSize)
                        error(message('bioinfo:dnds:NotValidWindowSize'))
                    end
                case 5 % ADJUSTSTOPS
                    adjustStops = bioinfoprivate.opttf(pval);
                    if isempty(verbose)
                        error(message('bioinfo:dnds:adjuststopsInputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
            end
        end
    end
end

if windowSize >= num_codons_original
    windowSize = num_codons_original;
    slidingWindow = false;
end

[translation, codons] = geneticCodeTable(gene_code);
N = nine_neighbors;

% gaps and ambiguous symbols are ignored, remove the whole codon
h = all(s1>=1,1) & all(s1<=4,1) & all(s2>=1,1) & all(s2<=4,1);

s1 = s1(:,h);
s2 = s2(:,h);
num_codons = sum(h);
num_nt = num_codons * 3;

seq1 = int2nt(s1);
seq2 = int2nt(s2);

% change nucleotide int triplets (1-4) to int codons (1-64)
s1 = (4.^(2:-1:0))*(double(s1)-1)+1;
s2 = (4.^(2:-1:0))*(double(s2)-1)+1;

if adjustStops
    %remove also any stop codons from the computations
    g = translation(s1)~='*' & translation(s2)~='*'; 
    s1 = s1(g);
    s2 = s2(g);
    seq1 = seq1(:,g);
    seq2 = seq2(:,g);
    h(h)=g;
    num_codons = sum(h);
    num_nt = num_codons * 3;
end

seq1 = seq1(:)';
seq2 = seq2(:)';

if num_codons<1
    error(message('bioinfo:dnds:SequenceTooShort'));
end

if verbose
    disp('DNDS: ')
    disp('Codons considered in the computations:')
    disp(seq1)
    disp(seq2)
    disp('Translations:')
    disp(sprintf('%c  ',nt2aa(seq1,'GeneticCode',gene_code)))
    disp(sprintf('%c  ',nt2aa(seq2,'GeneticCode',gene_code)))
end

if(strcmp(method , 'NG'))

    % Calculate the average potential syn and nonsyn sites for every codon
    % comparison in the sequences:
    s =  potential_synon_sites(N, translation);
    
    if adjustStops
        z = sum(N(:,translation~='*'),2);
        s = s * (num_codons*18/sum(z(s1)+z(s2)));
    end
    
    all_s = (s(s1(:))+s(s2(:)))/2;
    all_n = 3 - all_s;

    % Calculate the syn and nonsyn differences for every codon comparison
    % in the sequences:
    all_sd = zeros(num_codons,1);
    all_nd = zeros(num_codons,1);
    j = 1;
    for i = 1:3:num_nt-2
        [all_sd(j) all_nd(j)] = synonymous_diff(seq1(1,i:i+2),seq2(1,i:i+2), gene_code, adjustStops);
        j = j + 1;
    end

    % reinsert the codons with gaps and ambiguous symbols that had
    % previously been removed

    all_s=insertValue(all_s, h, 0);
    all_sd=insertValue(all_sd, h, 0);
    all_n=insertValue(all_n, h, 0);
    all_nd=insertValue(all_nd, h, 0);

    if slidingWindow
        win = ones(windowSize,1);

        S = convn(all_s,win,'valid');
        %N = (windowSize*3) - S;
        N = convn(all_n,win,'valid');

        Nd = convn(all_nd,win,'valid');
        Sd = convn(all_sd,win,'valid');

        S(S==0) = 1;        % rare case, it happens when Sd is also zero.
    else
        S = sum(all_s);
        N = num_nt - S;

        Sd = sum(all_sd);
        Nd = sum(all_nd);

        if S==0, S = 1; end % rare case, it happens when Sd is also zero.
    end

    pS = Sd ./ S;  % Nei&Gojobori (1)
    pN = Nd ./ N;  % Nei&Gojobori (2)

    pS(pS >= 3/4 ) = NaN;
    pN(pN >= 3/4 ) = NaN;

    if all(isnan(pS))
        % When using sliding windows, the warning is given only if all
        % positions are problematic.
        warning(message('bioinfo:dnds:psLarge'))
        if slidingWindow
            d_hat_s = nan(numel(pS),1);
            var_d_hat_s = nan(numel(pS),1);
        else
            d_hat_s = NaN;
            var_d_hat_s = NaN;
        end
    else
        d_hat_s = -(3/4)*log((1 - ((4/3)* pS))); % Nei&Gojobori (3)
        if slidingWindow
            var_p_hat_s=zeros(numel(pS),1);
            for w=1:numel(pS)
                var_p_hat_s(w) = sum(((all_sd(w:w+windowSize-1) - (all_s(w:w+windowSize-1)* pS(w))) .^2) / (S(w) .^2));
            end
        else
            var_p_hat_s = sum((( all_sd - (all_s * pS)).^2) / S^2); % Nei&Kumar (4.4)
        end
        var_d_hat_s = var_p_hat_s ./ ((1 - ((4/3)*pS)).^2); % Nei&Kumar (4.3)
    end

    if all(isnan(pN))
        % When using sliding windows, the warning is given only if all
        % positions are problematic.
        warning(message('bioinfo:dnds:pnLarge'))
        if slidingWindow
            d_hat_n = nan(numel(pN),1);
            var_d_hat_n = nan(numel(pN),1);
        else
            d_hat_n = NaN;
            var_d_hat_n = NaN;
        end
    else
        d_hat_n = -(3/4)*log((1 - ((4/3)* pN))); % Nei&Gojobori (3)
        if slidingWindow
            var_p_hat_n = zeros(numel(pN),1);
            for w=1:numel(pN)
                var_p_hat_n(w) = sum(((all_nd(w:w+windowSize-1) - (all_n(w:w+windowSize-1)* pN(w))) .^2) / (N(w) .^2)); %pf for sliding window
            end
        else
            var_p_hat_n = sum((( all_nd - (all_n * pN)).^2) / N^2); % Nei&Kumar (4.4)
        end
        var_d_hat_n = var_p_hat_n ./ ((1 - ((4/3)*pN)).^2);     % Nei&Kumar (4.3)
    end
    
    d_hat_s(d_hat_s==0) = 0; %remove signed zeros
    d_hat_n(d_hat_n==0) = 0; %remove signed zeros

else %LWL & PBL

    deg =  degeneracy(N, translation, codons);

    p_all=zeros(3,num_codons); %transition counts
    q_all=zeros(3,num_codons); %transversion counts
    l_all=zeros(3,num_codons); %average degeneracy

    for i = 1:3:num_nt-2
        [pqTemp lTemp] = trans_trav_deg_sub_count(seq1(1,i:i+2),seq2(1,i:i+2),deg(s1((i+2)/3),:),deg(s2((i+2)/3),:));
        p_all(1:3,(i+2)/3)=pqTemp(:,1); %p_all(x,y)=transition counts for position x in codon y
        q_all(1:3,(i+2)/3)=pqTemp(:,2); %q_all(x,y)=transversion counts for position x in codon y
        l_all(1:3,(i+2)/3)=lTemp;       % degeneracy
    end

    if slidingWindow

        %reinsert codons with ambiguous symbols and gaps that had
        %previously been removed
        t1=zeros(3,num_codons_original);
        t2=zeros(3,num_codons_original);
        t3=zeros(3,num_codons_original);

        for i = 1:3
            t1(i,:)=insertValue(p_all(i,:), h, 0);
            t2(i,:)=insertValue(q_all(i,:), h, 0);
            t3(i,:)=insertValue(l_all(i,:), h, 0);
        end
        p_all=t1;
        q_all=t2;
        l_all=t3;

        win = ones(1,windowSize);
        L = convn(l_all,win,'valid');
        L(L==0) = 1; % rare case, it happens when Q and P are also zero.

        P = convn(p_all,win,'valid') ./L;
        Q = convn(q_all,win,'valid') ./L;

    else
        L=sum(l_all,2);
        L(L==0) = 1; % rare case, it happens when Q and P are also zero.

        P=sum(p_all,2)./L;% proportion of transitional nucleotide differences at the ith class (i=0,2,4)
        Q=sum(q_all,2)./L;% proportion of transversional nucleotide differences at the ith class (i=0,2,4)
    end
    
    outLength = size(L,2);
    
    P((2*P)>=(1-Q)) = nan;
    Q(Q>=1/2) = nan;
    
    PLargeWarn = all(any(isnan(P),1));
    QLargeWarn = all(any(isnan(Q),1));
    % When using sliding windows, the warning is given only if all
    % positions are problematic. 
    if PLargeWarn || QLargeWarn
        if PLargeWarn
            warning(message('bioinfo:dnds:PLarge'))
        end
        if QLargeWarn
            warning(message('bioinfo:dnds:QLarge'))
        end
        if slidingWindow
           d_hat_n = nan(outLength,1);
           d_hat_s = nan(outLength,1);
           var_d_hat_n = nan(outLength,1);
           var_d_hat_s = nan(outLength,1); 
        else
           d_hat_n = NaN;
           d_hat_s = NaN;
           var_d_hat_n = NaN;
           var_d_hat_s = NaN; 
        end
        return
    end

    % Nei&Kumar (4.9)
    a = 1./(1-2*P-Q);
    b = 1./(1-2*Q);
    A = (1/2) .* log(a) - (1/4) .* log(b);
    B = (1/2) * log(b);
    c = (a - b)./2;
    d = b+c;
    varA = (((a.^2) .* P) + ((c.^2) .* Q) - ((a.*P) + (c.*Q)).^2)./L;
    varB = ((b.^2) .* Q .* (1 - Q))./L;
    varK = (((a.^2).*P) + ((d.^2).*Q) - (((a.*P) + (d.*Q)).^2)) ./L;

    if(strcmp(method, 'LWL'))
        d_hat_s = (3*((L(2,:).*A(2,:)) + (L(3,:) .* (A(3,:) + B(3,:))))) ./ (L(2,:) +(3* L(3,:)));
        d_hat_s(d_hat_s<0) = 0;
        d_hat_n = (3*((L(1,:) .* (A(1,:) + B(1,:))) + (L(2,:) .* B(2,:)))) ./ ((3*L(1,:)) + (2*L(2,:)));
        d_hat_n(d_hat_n<0) = 0;
        var_d_hat_s = (9*(((L(2,:).^2).*varA(2,:)) + ((L(3,:).^2).*varK(3,:))))./((L(2,:) + (3*L(3,:))).^2);
        var_d_hat_n = (9*(((L(2,:).^2).*varB(2,:)) + ((L(1,:).^2).*varK(1,:))))./(((2*L(2,:)) + (3*L(3,:))).^2);
    else % if(strcmp(method,'PBL'))
        d_hat_s = (((L(2,:) .* A(2,:)) + (L(3,:) .* A(3,:))) ./ (L(2,:) + L(3,:))) + B(3,:);
        d_hat_s(d_hat_s<0) = 0;
        d_hat_n = A(1,:) + (((L(1,:) .* B(1,:)) + (L(2,:) .* B(2,:))))./ (L(1,:) + L(2,:));
        d_hat_n(d_hat_n<0) = 0;
        var_d_hat_s = ((((L(2,:).^2).*varA(2,:))+((L(3,:).^2).*varA(3,:)))./((L(2,:)+L(3,:)).^2)+varB(3,:)) -...
            (((b(3,:).*Q(3,:)).*((2*a(3,:).*P(3,:))-c(3,:).*(1-Q(3,:))))./(L(2,:)+L(3,:)));
        var_d_hat_n = ((varA(1,:)+(((L(1,:).^2).*varB(1,:))+((L(2,:).^2).*varB(2,:))))./((L(1,:)+L(2,:)).^2)) - ...
            (((b(1,:).*Q(1,:)).*((2*a(1,:).*P(1,:))-(c(1,:).*(1-Q(1,:)))))./(L(1,:)+L(2,:)));
    end
    d_hat_s=d_hat_s';
    d_hat_n=d_hat_n';
    var_d_hat_s=var_d_hat_s';
    var_d_hat_n=var_d_hat_n';
    d_hat_s(d_hat_s==0) = 0; %remove signed zeros
    d_hat_n(d_hat_n==0) = 0; %remove signed zeros
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Common sub-functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------
function [translation, codons] = geneticCodeTable(gene_code)
codons = int2nt([ceil(1/16:1/16:4);repmat(ceil(1/4:1/4:4),1,4);repmat(1:4,1,16)])';
translation = nt2aa(codons,'GeneticCode',gene_code,'AlternativeStartCodons',false);

%-------------------------------------------
function N = nine_neighbors
% Initializes a 64x64 logical matrix indicating neighbors for each codon,
% neighbors are those that have only one nucleotide changed
v = eye(4,'uint8')==1;
u = ~v;
u = [u v v v; v u v v; v v u v; v v v u];
v = eye(16,'uint8')==1;
N = [u v v v; v u v v; v v u v; v v v u];

%-------------------------------------------
function [w]=insertValue(v,index,val)
% Insert value val into vector v at positions specified by logical vector
% index. Logical zeros indicate positions where the insertion should occur.

toInsert=find(index==0); % positions where val should be inserted
if(~isempty(toInsert))
    curr=toInsert(1);
    w=zeros(numel(v)+ numel(toInsert),1);
    w(1:toInsert(1)-1)=v(1:toInsert(1)-1);

    for idx=2:numel(toInsert)
        curr=toInsert(idx);
        prev=toInsert(idx-1);
        incr=curr-prev-1;
        if (incr> 0)
            w(prev+1:prev+incr)=v(prev+1-idx+1:prev+incr-idx+1);
        end
        w(prev)=val;
    end

    w(curr)=val;

    if(isempty(idx)) % only one codon to reinsert
        w(curr+1:end)=v(curr:end);
    else
        w(curr+1:end)=v(prev+incr-idx+1+1:end);
    end

    if(size(v,1)<size(v,2))
        w=w';
    end

else
    w=v;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NG sub-functions  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------
function syn  = potential_synon_sites(N,A)
% Finds the number potential synonymous changes for each possible
% codon, returns a 64x1 vector
Aint  = aa2int(A);
[i,j] = find(N);
N(N)  = Aint(i)==Aint(j);
syn   = sum(N,2)/3;

%------------------------------------------
function [sd nd] = synonymous_diff(codon1, codon2, gene_code, adjustStops)
% finds all possible paths from codon1 to codon2 and counts the number of
% synonymous changes and nonsynonymous changes. Note: no special case for
% nonsense mutations
nindx = find(codon1 ~= codon2);
sd = 0;
nd = 0;
if(~isempty(nindx))
    order = perms(nindx);
    [paths changes] = size(order);
    first = repmat(codon1,paths,1);
    second = repmat(codon2,paths,1);
    A = first;
    for i = 1:changes -1
        ind1 = sub2ind(size(first), (1:numel(order(:,1)))', order(:,i));
        first(ind1) = second(ind1);
        A = [A first]; %#ok<AGROW>
    end
    A = [A second];
    trsl = nt2aa(reshape(A',1,numel(A)),'GeneticCode',gene_code,'AlternativeStartCodons',false);
    B  = reshape(trsl,changes+1,paths)';
    if adjustStops
        B = B(~any(B=='*',2),:);
        paths = size(B,1);
    end
    syn = 0;
    non = 0;
    for i = 1:changes
        s = B(:,i) == B(:,i+1);
        syn = syn+sum(s);
        non = non+sum(~s);
    end
    sd = syn/paths;
    nd = non/paths;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LWL and PBL sub-functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------
function [PQ L] = trans_trav_deg_sub_count(codon1, codon2,deg1,deg2)
% tabulates the number of transitions and transversions at each site
% according to the codons degeneracy.
diffs = (codon1 ~= codon2);
PQ = zeros(5,2);
L =  zeros(3,1);
L(1) = L(1) + (sum([deg1 deg2] == 0) /2);
L(2) = L(2) + (sum([deg1 deg2] == 2) /2);
L(3) = L(3) + (sum([deg1 deg2] == 4) /2);
if(any(diffs))
    trans = is_transition(codon1,codon2);
    for i = 1:length(trans)
        if(trans(i) == 1)
            %transition
            PQ(deg1(i) + 1 ,1) = PQ(deg1(i) + 1 ,1) + .5;
            PQ(deg2(i) + 1 ,1) = PQ(deg2(i) + 1 ,1) + .5;
        elseif(trans(i) == 2)
            %transversion
            PQ(deg1(i) + 1 ,2) = PQ( deg1(i) + 1,2) + .5;
            PQ(deg2(i) + 1 ,2) = PQ( deg2(i) + 1,2) + .5;
        end
    end
end
PQ([2 4],:) = []; %reshape  PQ values

%------------------------------------------------------------
function  [deg] =  degeneracy(N, A, codons)
% gets degeneracy for each site returns its value as row vector of 3
% numbers {0,2,4}
% 4-fold degenerate all possible changes at each site is synonymous
% 2-fold degenerate 1 of 3 possible changes at each site is synonymous
% 0-fold degenerate all possible changes at each site is nonsynonymous
deg = zeros(64,3);
for codonindex = 1:64
    codon = codons(codonindex,:);
    neighs = N(codonindex,:);
    amino = A(codonindex);
    intcodons = codons(neighs,:);
    diffs = [codon(1)~= intcodons(:,1) codon(2) ~= intcodons(:,2) codon(3) ~= intcodons(:,3)];
    neightrans = A(neighs);
    stops = neightrans=='*';
    neightrans(stops) = ' ';
    index = amino==neightrans;
    degenerate = [sum(index&diffs(:,1)) sum(index&diffs(:,2)) sum(index&diffs(:,3))];
    onesthrees = (degenerate == 1 | degenerate == 3);
    degenerate(onesthrees) = degenerate(onesthrees) + 1;
    deg(codonindex,:) = degenerate;
end

%------------------------------------------------------------
function trans = is_transition(codon1, codon2)
% expects equal sized codon1 and codon2
% returns sequence of numerics
% 0 for equal bases
% 1 transition
% 2 transversion
A = codon1 ~= codon2;
tt =  (rem(nt2int(codon2(A)),2) ~= rem(nt2int(codon1(A)),2))+1;
trans = zeros(1, length(A));
trans(A) = tt;
