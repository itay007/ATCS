function out = cleavelookup(varargin)
%CLEAVELOOKUP displays cleavage rules of enzymes or compounds.
%
%   CLEAVELOOKUP displays a table of all abbreviation codes, cleavage
%   positions, cleavage patterns and full names of enzymes and compounds
%   for which cleavage rules are specified.
%
%   CLEAVELOOKUP('CODE', CODE) displays the cleavage position, cleavage
%   pattern and full name of the enzyme or compound specified by the
%   abbreviation CODE.
%
%   CLEAVELOOKUP('NAME', NAME) displays the cleavage position, cleavage
%   pattern and abbreviation code of the enzyme or compound whose full name
%   is specified by NAME. 
%
%   Examples:
%
%      cleavelookup('name', 'CASPASE 1')
%      cleavelookup('code', 'CASP1')
%
%   See also CLEAVE, REBASECUTS, RESTRICT.

%   Reference:
%   http://ca.expasy.org/tools/peptidecutter/peptidecutter_enzymes.html

%   Copyright 2008-2010 The MathWorks, Inc.



%=== check inputs
bioinfochecknargin(nargin,0,mfilename);
[code, name] = parse_inputs(varargin{:});

%=== read table with code, position, pattern and full name
out = '';
t = localClvTable;

%=== main
if nargin == 0  
    out = sprintf('%-9s%-9s%-70s%-31s\n','Code','Position','Pattern','Full Name');
     for i=1:length(t)
         %out = sprintf('%s%-9s%-7s%-52s%-31s\n', out, t{i,1}, num2str(t{i,2}), t{i,3}, t{i,4});
         out = sprintf('%s%-9s%-9s%-70s%-31s\n', out, t{i,1}, num2str(t{i,2}), t{i,3}, t{i,4});
     
     end
else
    if ~isempty(code)
        c = strmatch(upper(code), t(:,1), 'exact');
        if isempty(c)
            error(message('bioinfo:cleavelookup:UnknownCode', code));
        else
            out = sprintf('%s\t%s\t%s\t\n', num2str(t{c,2}), t{c,3}, t{c,4});
        end
    end
    if ~isempty(name)
        name = regexprep(name, '\s{2,}', ' '); % remove extra spaces
        c = strmatch(upper(name), t(:,4), 'exact');
        if isempty(c)
            error(message('bioinfo:cleavelookup:UnknownName', name));
        else
            out = sprintf('%s\t%s\t%s\t\n', num2str(t{c,2}), t{c,3}, t{c,1});
        end
    end
end
   

%--------------------------------------------------------------------------
function t = localClvTable
% Returns cell array with codes, positions, patterns and names.

t = cell(34,4);
t(1,:)  = {'ARG-C',  1, 'R',                        'ARG-C PROTEINASE'};
t(2,:)  = {'ASP-N',  2, 'D',                        'ASP-N ENDOPEPTIDASE'};
t(3,:)  = {'BNPS',   1, 'W',                        'BNPS-SKATOLE'};
t(4,:)  = {'CASP1',  1, '(?<=[FWYL]\w[HAT])D(?=[^PEDQKR])',  'CASPASE 1'};	
t(5,:)  = {'CASP2',  1, '(?<=DVA)D(?=[^PEDQKR])',            'CASPASE 2'};
t(6,:)  = {'CASP3',  1, '(?<=DMQ)D(?=[^PEDQKR])',            'CASPASE 3'};
t(7,:)  = {'CASP4',  1, '(?<=LEV)D(?=[^PEDQKR])',            'CASPASE 4'};
t(8,:)  = {'CASP5',  1, '(?<=[LW]EH)D',                      'CASPASE 5'};
t(9,:)  = {'CASP6',  1, '(?<=VE[HI])D(?=[^PEDQKR])',         'CASPASE 6'};
t(10,:) = {'CASP7',  1, '(?<=DEV)D(?=[^PEDQKR])',            'CASPASE 7'};
t(11,:) = {'CASP8',  1, '(?<=[IL]ET)D(?=[^PEDQKR])',         'CASPASE 8'};
t(12,:) = {'CASP9',  1, '(?<=LEH)D',                         'CASPASE 9'};
t(13,:) = {'CASP10', 1, '(?<=IEA)D',                         'CASPASE 10'};
t(14,:) = {'CH-HI',  1, '([FY](?=[^P]))|(W(?=[^MP]))',       'CHYMOTRYPSIN-HIGH SPECIFICITY'};
t(15,:) = {'CH-LO',  1, '([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))', 'CHYMOTRYPSIN-LOW SPECIFICITY'};
t(16,:) = {'CLOST',  1, 'R',                        'CLOSTRIPAIN'};
t(17,:) = {'CNBR',   1, 'M',                        'CNBR'};
t(18,:) = {'ENTKIN', 1, '(?<=[DN][DN][DN])K',       'ENTEROKINASE'};
t(19,:) = {'FACTXA', 1, '(?<=[AFGILTVM][DE]G)R',    'FACTOR XA'};
t(20,:) = {'FORMIC', 1, 'D',                        'FORMIC ACID'};
t(21,:) = {'GLUEND', 1, 'E',                        'GLUTAMYL ENDOPEPTIDASE'};
t(22,:) = {'GRANB',  1, '(?<=IEP)D',                'GRANZYME B'};
t(23,:) = {'HYDROX', 1, 'N(?=G)',                   'HYDROXYLAMINE'};
t(24,:) = {'IODOB',  1, 'W',                        'IODOSOBENZOIC ACID'};
t(25,:) = {'LYSC',   1, 'K',                        'LYSC'};
t(26,:) = {'NTCB',   1, 'C',                        'NTCB'};
t(27,:) = {'PEPS',   1, '((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|((?<=[^HKR][^P])[FLWY](?=\w[^P]))',   'PEPSIN PH = 1.3'};
t(28,:) = {'PEPS2',  1, '((?<=[^HKR][^P])[^R](?=[FL][^P]))|((?<=[^HKR][^P])[FL](?=\w[^P]))',  'PEPSIN PH > 2'};
t(29,:) = {'PROEND', 1, '(?<=[HKR])P(?=[^P])',      'PROLINE ENDOPEPTIDASE'};
t(30,:) = {'PROTK',  1, '[AEFILTVWY]',              'PROTEINASE K'};
t(31,:) = {'STAPHP', 1, '(?<=[^E])E',               'STAPHYLOCOCCAL PEPTIDASE I'};
t(32,:) = {'THERMO', 1, '[^DE](?=[AFILMV])',        'THERMOLYSIN'};
t(33,:) = {'THROMB', 1, '((?<=\w\wG)R(?=G))|((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))', 'THROMBIN'};
t(34,:) = {'TRYPS',  1, '((?<=\w)[KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))', 'TRYPSIN'};

%--------------------------------------------------------------------------
function [code, name] = parse_inputs(varargin)
% Parse input PV pairs.
     
%=== Check for the right number of inputs
if rem(nargin,2)== 1
    error(message('bioinfo:cleavelookup:IncorrectNumberOfArguments', mfilename));
end

%=== Defaults
code = '';
name = '';

%=== Allowed inputs
okargs = {'code','name'};

for j=1:2:nargin
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1 % code
            if ~isempty(pval) && ischar(pval)
                code = pval;
            else
                error(message('bioinfo:cleavelookup:InvalidCode'));
            end
        case 2 % name
            if  ~isempty(pval) && ischar(pval)
                name = pval;
            else
                error(message('bioinfo:cleavelookup:InvalidName'));
            end
    end      
end
