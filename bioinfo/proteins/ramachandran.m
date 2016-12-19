function out =  ramachandran(pdb_struct, varargin)
%RAMACHANDRAN computes torsion angles and generates Ramachandran plot.
%
%   RAMACHANDRAN(PDBID) generates the Ramachandran plot for the protein
%   structure specified by PDB code PDBID.
%
%   RAMACHANDRAN(PDBFILE) generates the Ramachandran plot for the protein
%   structure represented in PDB file PDBFILE. 
%
%   RAMACHANDRAN(PDBSTRUCT) generates the Ramachandran plot for the protein
%   structure stored in PDBSTRUCT, where PDBSTRUCT is a MATLAB structure
%   obtained by using PDBREAD or GETPDB.   
%
%   RAMACHANDRAN generates a plot of the torsion angles PHI (torsion angles
%   between the 'C-N-CA-C' atoms) and the torsion angles PSI (torsion
%   angles between the 'N-CA-C-N' atoms) of a protein.
%
%   RAMASTRUCT = RAMACHANDRAN(...) returns a MATLAB structure or array of
%   structures with the following fields:
%      Angles - torsion angles PHI, PSI and OMEGA. 
%      ResidueNum - residue sequence numbers. 
%      ResidueName - residue names. 
%      Chain - chain ID in the PDB structure.
%      HPoints - handle to the data points in the plot.
%
%   The torsion angles are in the same order as the residue numbers and the
%   residue names. The number of elements in RAMASTRUCT is equal to the
%   number of chains considered. 
% 
%   RAMACHANDRAN(..., 'CHAIN', CHAINID) computes the torsion angles and
%   generates the Ramachandran plot for residues in chain(s) CHAINID. Valid
%   choices are: 'ALL' to indicate all available chains, a string
%   specifying the chain ID to consider, or a cell array of strings
%   specifying the list of chain IDs to consider. CHAINID is case
%   sensitive. Default is 'ALL'.
% 
%   RAMACHANDRAN(..., 'PLOT', PLOTVALUE) specifies how the torsion angles
%   should be plotted with respect to their chains. Valid choices are: 
%      'none' - do not generate any plot.
%      'separate' - plot each chain separately. 
%      'combined' - plot all considered chains in a single plot. 
%   Default is 'combined'.
% 
%   RAMACHANDRAN(..., 'MODEL', MODELNUM) specifies the PDB structural
%   model to consider. Default is 1.
%   
%   RAMACHANDRAN(..., 'GLYCINE', TF) controls the highlighting of glycine
%   residues with a circle in the plot. Default is false.
%
%   RAMACHANDRAN(..., 'REGIONS', TF) controls the drawing of reference
%   regions in the plot. The default reference regions are: core
%   right-handed alpha, core beta, core left-handed alpha, allowed. The
%   boundaries of default reference regions are based on the calculations
%   by Morris et al., 1992. Default is false.
%
%   RAMACHANDRAN(..., 'REGIONDEF', RDEF) specifies the names, colors and
%   boundaries of the custom reference regions in the Ramachandran plot(s).
%   RDEF is a MATLAB structure containing the following fields:
%      Name - string specifying a name for a region.
%      Color - string or RGB vector specifying a color for a region.
%      Patch - 2xN matrix of values defining the boundaries of a region. N
%      is the number of data points needed to define the region. 
%
%   Examples:
% 
%   % Generate the Ramachandran plot for the human serum albumin
%   % complexed with octadecanoic acid.
%   ramachandran('1E7I')
%
%   % Compute the torsion angles and the Ramachandran plot for chain A of 
%   % the human growth hormone, represented in the PDB file 1a22.pdb. 
%   getpdb('1a22', 'tofile', '1a22.pdb');
%   a = ramachandran('1a22.pdb', 'chain', 'A')
%
%   % Generate a separate Ramachandran plot for each chain in the pdb file 
%   1a22.pdb, highlighting the Ramachandran regions and the glycine
%   residues. 
%   ramachandran('1a22.pdb', 'plot', 'separate', 'glycine', true, 'regions', true)
% 
%   See also GETPDB, MOLVIEWER, PDBDISTPLOT, PDBREAD, PROTEINPROPPLOT.

%   References:
%   A. Morris, M. MacArthur, G. Hutchinson, and J. Thornton. Stereochemical
%   Quality of Protein Structure Coordinates. Proteins: Structure,
%   functions, and Genetics (1992) 12:345-364

% Copyright 2003-2010 The MathWorks, Inc.



%=== check inputs
bioinfochecknargin(nargin,1,mfilename)
[chain, plotType, model, regionFlag, glycineFlag, heteroFlag, regData] = ...
    parse_inputs(varargin{:});

%=== get pdb structure
try
    if (~isstruct(pdb_struct))
        if exist(pdb_struct,'file') % read it from the file
            pdb_struct = pdbread(pdb_struct);
        else % get it from PDB
            pdb_struct = getpdb(pdb_struct);
        end
    else % get it from structure
        pdb_struct = convertpdbstruct(pdb_struct, mfilename);
    end
catch theErr
    if ~isempty(strfind(theErr.identifier,'getpdb')) 
        rethrow(theErr);
    end
    error(message('bioinfo:ramachandran:IllegalInput'));
end
 
%=== Get the model from pdbstruct
if isfield(pdb_struct, 'Model')
    try
        model_struct = pdb_struct.Model(model);
    catch allExceptions
        error(message('bioinfo:ramachandran:InvalidModel'));
    end
else
    error(message('bioinfo:ramachandran:NotAtomicCoordinates'));
end

%=== Get hetero atom records if wanted
if heteroFlag && isfield(model_struct, 'HeterogenAtom')
    hetAtoms = model_struct.HeterogenAtom;
else
    hetAtoms = [];
end

%=== Get atom records
if isfield(model_struct, 'Atom')
    atoms = model_struct.Atom;
else
    atoms = [];
end

%=== Obtaining the serial numbers for the Atoms and Heterogen Atoms
a = zeros(1,length(atoms)+length(hetAtoms));
for i = 1:length(atoms)
    a(i) = atoms(i).AtomSerNo;
end
for j = 1:length(hetAtoms)
    a(i+j) = hetAtoms(j).AtomSerNo;
end

pdb_struct.NewAtom = [atoms hetAtoms];

[sorted_SerNo, sorted_index] = sort(a); 
pdb_struct.NewAtom = pdb_struct.NewAtom(sorted_index);
clear sorted_SerNo sorted_index a


%=== Extract atoms needed to compute torsion angles
[atom_name, atom_resSeq, atom_coord, atom_chain, atom_resName] = localGetAtoms(pdb_struct);


%=== Select chain(s) to consider
[chains, chainStart, chainStop] = localGetChains(atom_chain, chain);
nChains = numel(chains);


%=== Compute torsion angles
[out, m] = localGetTorsionAngles(atom_name, atom_resSeq, atom_resName, ...
    atom_coord, chains, chainStart, chainStop);


%=== Plot if required
if plotType > 1
    
    if ismac 
        glyMarkerSize = 8;
    else 
        glyMarkerSize = 6;
    end
    
    if nChains > 1 % multi chain
        if plotType == 3 % combine into array the data to be plotted
            M = max(m); % maximum number of residue among all chains
           
            phi  = ones(M, nChains) * NaN;
            psi  = ones(M, nChains) * NaN;
            resn = ones(M, nChains) * NaN;
            res  = cell(M, nChains);
            
            for c = 1:nChains
                phi(1:m(c),c) = out(c).Angles(:,1);
                psi(1:m(c),c) = out(c).Angles(:,2);
                res(1:m(c),c) = out(c).ResidueName;
                resn(1:m(c),c) = out(c).ResidueNum;
            end
            
            %=== plot multiple chains, combined
            hFig = figure();
            set(hFig, 'tag', 'Ramachandran');
            ts = pdb_struct.Header.idCode; % title
            
            if regionFlag
                localPlotRamaRegions(regData, hFig);
            end
            
            line = localPlotRama(phi, psi, res, resn, chains, ts, hFig);
            
            for c = 1:nChains
                out(c).HPoints = line(c);
            end
             
            title(get(out(c).HPoints, 'parent'), ts);
            
            %=== highlight glycine
            if glycineFlag
                figure(hFig); 
                hold on;
                for c = 1:nChains
                    g = find(strcmp(out(c).ResidueName, 'GLY'));
                    if ~isempty(g)
                        hgly = plot(phi(g,c), psi(g,c), 'ko', ...
                            'MarkerSize', glyMarkerSize);
                        set(hgly,'HitTest','off', 'DisplayName', 'Glycines');
                    end
                end
                hold off;
            end
            
        else
            %=== plot multiple chains in separate figures
            for c = 1:nChains
                             
                hFig = figure();
                set(hFig, 'tag', 'Ramachandran');
                               
                ts = [pdb_struct.Header.idCode ' Chain ' chains{c}]; %title
                if regionFlag
                    localPlotRamaRegions(regData, hFig);
                end
                
                phi = out(c).Angles(:,1);
                psi = out(c).Angles(:,2);
                
                out(c).HPoints = localPlotRama(phi, psi, out(c).ResidueName, ...
                    out(c).ResidueNum, chains{c}, ts, hFig);
                
                title(get(out(c).HPoints, 'parent'), ts); 
                
                %=== highlight glycines
                if glycineFlag
                    g = find(strcmp(out(c).ResidueName, 'GLY'));
                    if ~isempty(g)
                        figure(hFig); 
                        hold on;
                        hgly = plot(phi(g), psi(g), 'ko', ...
                            'MarkerSize', glyMarkerSize, 'DisplayName', 'Glycines');
                        set(hgly, 'HitTest','off')
                        hold off;

                    end
                end

            end
            
        end

    else %=== plot single chain
        hFig = figure();
        set(hFig, 'tag', 'Ramachandran');
        ts = [pdb_struct.Header.idCode ' Chain ' chains{1}]; %title

        if regionFlag
            localPlotRamaRegions(regData, hFig);
        end

        phi = out.Angles(:,1);
        psi = out.Angles(:,2);
        out.HPoints = localPlotRama(phi, psi, out.ResidueName, ...
            out.ResidueNum, chains{1}, ts, hFig);
        
        if glycineFlag
            g = find(strcmp(out.ResidueName, 'GLY'));
            if ~isempty(g)
                figure(hFig); 
                hold on;
                hgly = plot(phi(g), psi(g), 'ko', ...
                    'MarkerSize', glyMarkerSize, 'DisplayName', 'Glycines');
                set(hgly, 'HitTest','off')
                hold off;
            end
        end
        title(get(out.HPoints, 'parent'), ts);
    end
    
end
    


%==========================================================================
% SUBFUNCTIONS
%==========================================================================

function HPoints = localPlotRama(phi, psi, res, resn, c, ts, hFig)
% Plot torsion angles.

figure(hFig);
hold on;
hAxis = gca;

%=== figure settings
axis([-180 180 -180 180])
axis square; box on;
hl = line([0 -180; 0 180], [-180 0; 180 0], 'Color', 'k', 'LineStyle', ':');
set(hl, 'HitTest', 'off'); % exclude line from displying cursor

%=== exclude two perpendicular lines from legend
set(get(get(hl(1),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off');
set(get(get(hl(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); 

lim = linspace(-180, 180, 9);
set(hAxis, 'Xtick', lim, 'YTick', lim);

set(hFig, 'NumberTitle', 'off', 'Name', ['Ramachandran Plot: ' ts]);
set(hAxis,'Layer', 'top')
xlabel('Phi (Degrees)');
ylabel('Psi (Degrees)');

%=== plot data
HPoints = plot(phi, psi, 'k.');
for i=1:numel(HPoints) % exclude from legend
    set(get(get(HPoints(i), 'Annotation'), 'legendInformation'), ...
    'IconDisplayStyle', 'off');
end
hold off;

%=== data cursor
datacursormode on;
dcm_obj = datacursormode(hFig);
set(dcm_obj,'UpdateFcn',{@myupdatefcn, phi, psi, res, resn, c})

%--------------------------------------------------------------------------

function datacursorLabel = myupdatefcn(obj, event_obj, phi, psi, res, resn, c) %#ok<INUSL>
% Display the data cursor.
 
k = 1; % default: phi and psi are vectors

tg = get(event_obj, 'Target');
if ~strcmpi(get(tg, 'Type'), 'patch')
    j = get(event_obj, 'DataIndex');
    if size(phi, 2) > 1
        pos = get(event_obj, 'Position');
        k = find(phi(j,:) == pos(1));
        c = c{k};
    end
    
    r = res{j,k};
    n = resn(j,k);
   
    datacursorLabel = {[r ' ' num2str(n)], ...
        ['Phi: ',num2str(phi(j,k),4)], ...
        ['Psi: ',num2str(psi(j,k),4)], ...
        ['Chain: ' c]};
else
        region = get(tg, 'DisplayName');
        datacursorLabel = {region};
end

%--------------------------------------------------------------------------
function localPlotRamaRegions(regData, hFig)
% Plot patches corresponding to reference Ramachandran regions as described in regData.

hp = zeros(1,numel(regData));
figure(hFig); 

for i = numel(regData):-1:1 

    % print only contours
     hp(i) = patch(regData(i).Patch(1,:), regData(i).Patch(2,:),...
         regData(i).Color, 'EdgeColor', regData(i).Color,...
        'DisplayName', regData(i).Name, 'FaceColor', 'none');
    
%    % print solid patches
%     hp(i) = patch(regData(i).Patch(1,:), regData(i).Patch(2,:), regData(i).Color, ...
%         'DisplayName', regData(i).Name, 'EdgeColor', 'none');
   
%     set(get(get(hp, 'Annotation'), 'legendInformation'),
%     'IconDisplayStyle', 'off')
end

% %== group
% if ischar(regData(1).Color) % string color assignments
%     rType = unique({regData.Color}); % distinct region types
% else
%    rType = unique(vertcat(regData.Color), 'rows'); % distinct region types
% end
% for  i = length(rType):-1:1
%     hg = hggroup;
%     f = strcmp(rType(i), {regData.Color});
%     set(hp(f), 'Parent', hg);
%     set(get(get(hg,'Annotation'),'LegendInformation'),...
%         'IconDisplayStyle','on'); % Include this hggroup in the legend
%     set(hg, 'DisplayName', ['Region Type ' num2str(i)]);
% end

%---------------------------------------------------------------------------

function [atom_name, atom_resSeq, atom_coord, atom_chain, atom_resName] = ...
    localGetAtoms(pdb_struct)
% Extract atoms 'C', 'CA', or 'N' and their coordinates (only these atoms
% are responsible for the torsion angles).

numAtoms = length(pdb_struct.NewAtom);
atom_name = cell(1,numAtoms); % name of atom (e.g., CA, N, etc.)
atom_resSeq = zeros(1,numAtoms); % residue number in the chain
atom_coord = zeros(numAtoms,3); % atom coordinates
atom_chain = cell(1,numAtoms); % chain
atom_resName = cell(1,numAtoms); % residue name (three-letter symbol)

k = 1;
for i = 1:numAtoms
    if (isequal(pdb_struct.NewAtom(i).AtomName, 'C') ||...
        isequal(pdb_struct.NewAtom(i).AtomName, 'N') ||...
        isequal(pdb_struct.NewAtom(i).AtomName, 'CA'))
        atom_coord(k,:) = [pdb_struct.NewAtom(i).X...
                           pdb_struct.NewAtom(i).Y...
                           pdb_struct.NewAtom(i).Z];
        atom_name{k} = pdb_struct.NewAtom(i).AtomName;
        atom_resSeq(k) = pdb_struct.NewAtom(i).resSeq;
        atom_chain{k} = pdb_struct.NewAtom(i).chainID;
        atom_resName{k} = pdb_struct.NewAtom(i).resName;
        k = k+1;
    end
end
atom_name(k:end) = [];
atom_resSeq(k:end) = [];
atom_coord(k:end,:) = [];
atom_chain(k:end) = [];
atom_resName(k:end) = [];


if isempty(atom_name)
    error(message('bioinfo:ramachandran:NotProteinStructure'));
end

%---------------------------------------------------------------------------

function [chains, chainStart, chainStop] = localGetChains(atom_chain, chain)
% Determine the chains to consider and the first and last atom for each
% chain.

[allChains, chainStart] = unique(atom_chain, 'first'); % first atom number in each chain
[allChains, chainStop] = unique(atom_chain, 'last'); % last atom number in each chain

n = numel(allChains);
if isempty(chain) 
    if n > 1
        warning(message('bioinfo:ramachandran:MultipleChains'));
        ch = (1:n)';
    else
        ch = 1;
    end
elseif strcmpi(chain, 'ALL')
    ch = (1:n)';
else
    w = warning;
    warning('off','bioinfo:seqmatch:StringNotFound'); % do not warn if no chain match (handled later on)
    ch = seqmatch(chain, allChains); % order of selected chains in the pdb
    warning(w);% restore warning state
end

if any(ch == 0) 
    error(message('bioinfo:ramachandran:InvalidChain')); 
else
    chains = {allChains{ch}};
    chainStart = chainStart(ch);
    chainStop = chainStop(ch);
end


%---------------------------------------------------------------------------

function [out, m] = localGetTorsionAngles(atom_name, atom_resSeq,atom_resName, ...
        atom_coord, chains, chainStart, chainStop)
% Create the output structure with torsion angles, name of residues,
% position of residue and hadle to the datapoints.

dummy.Angles = [];
dummy.ResidueNum = [];
dummy.ResidueName = '';
dummy.Chain = '';
dummy.HPoints = [];
out = repmat(dummy, 1, numel(chains));
m = zeros(1,numel(chains)); % vector holding the max residue number in each chain

for c = 1:numel(chains)
    
    i = chainStart(c); % atom index

    %=== Preprocessing the Data
    % Checking to see if the initial pattern is N-CA-C-N. In this case we
    % reject this data point as this would mean that both PHI and PSI are
    % not sharing the same CA (PHI would not be defined for the first
    % residue, while PSI would be).

    if ((isequal(atom_name{i}, 'N')) && (isequal(atom_name{i+1}, 'CA')))
        i = i + 2;
    end

    %=== get residue numbers
    min_atom_resSeq = min(atom_resSeq(chainStart(c):chainStop(c)));
    max_atom_resSeq = max(atom_resSeq(chainStart(c):chainStop(c)));
    offset = 0;
    if min_atom_resSeq <= 0
        offset = -(min_atom_resSeq) + 1;
    end
    out(c).ResidueNum = (1-offset:max_atom_resSeq)';
    m(c) = length(out(c).ResidueNum);

    %=== compute torsion angles
    %numResidues = max(atom_resSeq) + offset;
    numResidues = max_atom_resSeq + offset;
    phi = nan(numResidues,1);
    psi = phi;
    omega = phi;
    
    %=== angles not defined unless 4 consecutive atoms in the chain
    while(i <= chainStop(c) - 3) 
        switch(atom_name{i})
            case 'C'
                % does the pattern satisfies for the phi torsion angle?
                if (isequal(atom_name{i+1}, 'N') &&...
                        (isequal(atom_name{i+2}, 'CA')) &&...
                        (isequal(atom_name{i+3}, 'C')))
                    % Are atoms are from adjacent amino acids? If not then
                    % torsion angle does not exist (is NaN)
                    if ((atom_resSeq(i) == atom_resSeq(i+1)-1) && ...
                            (atom_resSeq(i+1) == atom_resSeq(i+2)) &&...
                            (atom_resSeq(i+1)== atom_resSeq(i+3)))
                        phi(atom_resSeq(i+1)+offset) = localCalculateTorsionAngle(atom_coord(i:i+3,:));
                    end

                end
            case  'N'
                % does the pattern satisfies for the psi torsion angle?
                if (isequal(atom_name{i+1}, 'CA') &&...
                        (isequal(atom_name{i+2}, 'C')) &&...
                        (isequal(atom_name{i+3}, 'N')))
                    if ((atom_resSeq(i) == atom_resSeq(i+1)) &&...
                            (atom_resSeq(i) == atom_resSeq(i+2)) &&...
                            (atom_resSeq(i) == atom_resSeq(i+3)-1))
                        psi(atom_resSeq(i)+offset) = localCalculateTorsionAngle(atom_coord(i:i+3,:));
                    end
                end
            case 'CA'
                %if (nargout > 0) % compute omega only if user asks for output
                    % does the pattern satisfies for the omega torsion angle?
                    if (isequal(atom_name{i+1}, 'C') &&...
                            (isequal(atom_name{i+2}, 'N')) &&...
                            (isequal(atom_name{i+3}, 'CA')))
                        if ((atom_resSeq(i) == atom_resSeq(i+1)) &&...
                                (atom_resSeq(i) == atom_resSeq(i+2)-1) &&...
                                (atom_resSeq(i) == atom_resSeq(i+3)-1))
                            omega(atom_resSeq(i)+offset) = localCalculateTorsionAngle(atom_coord(i:i+3, :));
                        end
                    end
               %end
        end
        i = i + 1;
    end
     
     %=== Converting the angles from Radians to Degrees
     phi = phi*180/pi;
     psi = psi*180/pi;
     out(c).Angles = [phi, psi, omega];

     %=== Get residue names
     resName = cell(numResidues, 1);
     un_resSeq = unique(atom_resSeq(chainStart(c):chainStop(c)));
     for i = 1:length(un_resSeq)
         idx = find(atom_resSeq==un_resSeq(i), 1);
         resName{un_resSeq(i)+ offset} = atom_resName{idx};
     end

     out(c).ResidueName = resName;
     out(c).Chain = chains{c};
end

%---------------------------------------------------------------------------

function [torsionAngle] = localCalculateTorsionAngle(coord)
% Evaluate the torsion angles.

p1 = coord(1,:);
p2 = coord(2,:);
p3 = coord(3,:);
p4 = coord(4,:);

a = p2 - p1;
b = p3 - p2;
c = p4 - p3;

a = a/norm(a,2);
b = b/norm(b,2);
c = c/norm(c,2);

b_cross_c = [b(2).*c(3) - b(3).*c(2);
    b(3).*c(1) - b(1).*c(3);
    b(1).*c(2) - b(2).*c(1)];

x = -sum(conj(a).*c) + ((sum(conj(a).*b))*(sum(conj(b).*c)));
y = sum(conj(a).*b_cross_c');
torsionAngle = localNewAngle(x,y);

%---------------------------------------------------------------------------

function ang = localNewAngle(x,y)
% Calculate the angle represented by (x,y). The angle lies between -pi and
% pi.

ang = 0; % This is the default value. In case y == 0 and x ~= 0, ang = 0.
if (x ~= 0) && (y~= 0)
    c = x./sqrt(x.^2 + y.^2);
    ang = sign(y)*acos(c);
elseif (x == 0)
    if (y > 0)
        ang = pi/2;
    elseif (y < 0)
        ang = -pi/2;
    end
end

%---------------------------------------------------------------------------


function [chain, plotType, model, regionFlag,  glycineFlag, heteroFlag, regData] ...
         = parse_inputs(varargin)
% Parse input PV pairs.
     
%=== Check for the right number of inputs
if rem(nargin,2)== 1
    error(message('bioinfo:ramachandran:IncorrectNumberOfArguments', mfilename));
end

%=== Allowed inputs
okargs = {'chain','plot','model','regions','glycine','hetero','regiondef'};

%=== Defaults
chain = '';          % chain to consider
model = 1;           % model to consider
regionFlag = false;  % display regions in plot
glycineFlag = false; % highlight Glycines
heteroFlag = false;  % consider hetero atoms
regData = feval('ramachandranRegions'); % default reference region information
plotType = 3;        % type of plot (none=1, separate=2 or combined=3)

for j=1:2:nargin
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1 % chain
            if iscell(pval)
                chain = pval;
            elseif ischar(pval)
                chain = {pval};
            else
                error(message('bioinfo:ramachandran:InvalidChainOption'));
            end
        case 2 % plot
            okplots = {'none', 'separate','combined'};
            plotType = find(strncmpi(pval, okplots, numel(pval)), 1);
            if isempty(plotType)
                error(message('bioinfo:ramachandran:InvalidPlotOption'));
            end
        case 3 % model
            if isnumeric(pval) && isscalar(pval)
                model = pval;
            else
                error(message('bioinfo:ramachandran:InvalidModelOption'));
            end
        case 4 % regions
             regionFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 5 % glycine
            glycineFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 6 % hetero
            heteroFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 7  % regiondef
            if isstruct(pval)
                regData = pval;
                regionFlag = localCheckRegionDef(regData);
            else
                error(message('bioinfo:ramachandran:InvalidRefRegion'));   
            end
    end
end

%--------------------------------------------------------------------------
function regionFlag = localCheckRegionDef(regData)
% make sure the structure passed with option regiondef is valid

if isfield(regData, 'Name') && isfield(regData, 'Color') && isfield(regData, 'Patch')
    regionFlag = true;
    charFlag = 0; % =1 if at least one color is given by a string label
    vecFlag = 0;  % =1 is at least one color is given as RGB vector

    %=== check for consistent color assignments
    for r = 1:numel(regData)
        charFlag = charFlag | ischar(regData(r).Color);
        vecFlag = vecFlag | isnumeric(regData(r).Color);
    end

    if charFlag && vecFlag
        error(message('bioinfo:ramachandran:InconsistentColorDef'));
    end

    if ~charFlag && ~vecFlag
        error(message('bioinfo:ramachandran:InvalidColorDef'));
    end

else
    error(message('bioinfo:ramachandran:MissingFieldRefRegion'));
end
