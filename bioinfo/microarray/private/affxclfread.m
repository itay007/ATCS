function out = affxclfread(filename)
%AFFXCLFREAD reads Affymetrix Chip Layout Files (CLF).
%
%   clfStruct = AFFXCLFREAD(FILE) reads an Affymetrix CLF file FILE
%   and creates a structure clfStruct. CLF files contain information about
%   the names chip layout.
%
%   CLF files have an optional 'Sequential' flag. If this is set then the
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


% The spec for CLF files is given here:
% https://www.affymetrix.com/support/developer/fusion/File_Format_CLF_aptv161.pdf

% open the file
fid = fopen(filename,'r');
c = onCleanup(@()fclose(fid));

if fid == -1
    error(message('bioinfo:affxclfread:badCLFFile', filename));
end

out = struct('ChipType','','LibSetName','','LibSetVersion','','CreateDate','',...
    'GUID','','ClfFormatVersion','','Rows',NaN,'Cols',NaN,'StartID',0,'EndID',0,'Order','',...
    'DataColNames','','Data',[]);


% The sequential field is optional. If it exists we can ignore the first
% column of the data.
isSequential = false;

% read the header
try
    pos = ftell(fid);
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
                case 'clf_format_version'
                    out.ClfFormatVersion = value;
                case 'rows'
                    out.Rows = str2double(value);
                case 'cols'
                    out.Cols = str2double(value);
                case 'sequential'
                    out.StartID = str2double(value);
                    isSequential = true;
                case 'order'
                    out.Order = value;
                case 'header0'
                    colNames = textscan(value,'%s','CollectOutput',true);
                    out.DataColNames = colNames{1}';
            end
        end
        pos = ftell(fid);
        headerLine = fgetl(fid);
    end
    fseek(fid,pos,-1);
    % now the data
    if isSequential % we don't need to store the data, just a way to calculate it
        funcStr = sprintf(...
            'out.Data = @(ProbeID)affxCLFSequentialData(ProbeID,%d,%d,%d,''%s'');',...
            out.Rows, out.Cols, out.StartID,out.Order);
        eval(funcStr);
        out.EndID = out.Rows*out.Cols+(out.StartID-1);
    else
        data = textscan(fid,'%d%d%d','delimiter','\t','CollectOutput',true);
        out.Data = data{1};
        out.StartID = double(min(out.Data(:,1)));
        out.EndID = double(max(out.Data(:,1)));
    end
    if isempty(out.ChipType)
         warning(message('bioinfo:affxclfread:CLFNoChipType', filename))
    end
catch allExceptions
    error(message('bioinfo:affxclfread:CannotReadCLFFile', filename))
end

function out = affxCLFSequentialData(probeID,rows,cols,start,order) %#ok<DEFNU>
% Function used to calculate the x and y position for a probeID when the
% file is deterministic. I.e. When the "sequential" field is present in the
% header.
probeID = probeID-start;

if strcmpi(order,'col_major')
    x = floor(probeID/rows);
    y = rem(probeID,rows);
else
    x = floor(probeID/cols);
    y = rem(probeID,cols);
end

out = [x(:),y(:)];
