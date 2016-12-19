function out = affxbgpread(filename)
%AFFXBGPREAD reads Affymetrix Chip Layout Files (BGP).
%
%   bgpStruct = AFFXBGPREAD(FILE) reads an Affymetrix BGP file FILE
%   and creates a structure bgpStruct. BGP files contain information about
%   the names chip layout.
%
%   BGP files have an optional 'Sequential' flag. If this is set then the
%   probe_id values are sequential in the file so there is no need to store
%   the values. Instead a function handle is returned that will calculate
%   the x and y probe position given a probe ID.
%
%   See also AFFYREAD, AGFEREAD, CELINTENSITYREAD, GPRREAD,
%   PROBELIBRARYINFO, PROBESETLINK, PROBESETLOOKUP, PROBESETPLOT,
%   PROBESETVALUES, SPTREAD.
%
%   GeneChip and Affymetrix are registered trademarks of Affymetrix, Inc.

% Copyright 2008-2009 The MathWorks, Inc.


% The spec for BGP files is given here:
% https://www.affymetrix.com/support/developer/fusion/File_Format_BGP_aptv161.pdf

% open the file
fid = fopen(filename,'r');
c = onCleanup(@()fclose(fid));

if fid == -1
    error(message('bioinfo:affxbgpread:badBGPFile', filename));
end

dataStruct = struct('probe_id',[],'probeset_id',[],	'type','','gc_count',[],...
    'probe_length',[],'interrogation_position',[],'probe_sequence','','atom_id',[],'x',[],'y',[]);

out = struct('ChipType','','LibSetName','','LibSetVersion','','CreateDate','',...
    'GUID','','ExecGUID','','ExecVersion','','Cmd','','Data',dataStruct);

delimiter = '\t';

% read the header
try
    headerLine = fgetl(fid);
    while(strncmpi(headerLine,'#',1))
        if(strncmpi(headerLine,'#%',2)) && length(headerLine)>=3 % ignore comments
            % could use regexp with split option but there is risk if we
            % get more than one = on a line. Go with safer strtok.
            [headerType, value] = strtok(headerLine(3:end),'=');
            value = value(2:end); % strip out the '='
            switch(headerType)
                case 'chip_type'
                    if isempty(out.ChipType)
                        out.ChipType = value;
                    else
                        if ischar(out.ChipType)
                            out.ChipType = {out.ChipType};
                        end
                        out.ChipType{end+1} = value;
                    end
                case 'lib_set_name'
                    out.LibSetName = value;
                case 'lib_set_version'
                    out.LibSetVersion = value;
                case 'create_date'
                    out.CreateDate = value;
                case 'guid'
                    out.GUID = value;
                case 'exec_version'
                    out.ExecVersion = value;
                case 'exec_guid'
                    out.ExecGUID = value;
                case 'cmd'
                    out.Cmd = value;
            end
        end
        headerLine = fgetl(fid);
    end
    
    % First read the column headers
    colHeaders = strread(headerLine,'%s','delimiter',delimiter);
    isNumericCol = true(numel(colHeaders),1);
    % now the data
    data = bioinfoprivate.bioReadMixedData(fid,isNumericCol,delimiter);
    
    %   data = textscan(fid,'%d%d%d','delimiter','\t','CollectOutput',true);
    for count = 1:numel(colHeaders)
        out.Data.(genvarname(colHeaders{count})) = data{count};
    end
    
    if isempty(out.ChipType)
        warning(message('bioinfo:affxbgpread:BGPNoChipType', filename))
    end
catch allExceptions
    error(message('bioinfo:affxbgpread:CannotReadBGPFile', filename))
end


