function aboutstr = bioinfoabout(varargin)%#ok
%BIOINFOABOUT About the Bioinformatics Toolbox.

tlbx = ver('bioinfo');
tlbx = tlbx(1);
aboutstr = sprintf('%s\nVersion %s\nCopyright 2003-%s, The MathWorks, Inc.', ...
    tlbx.Name, tlbx.Version, datestr(tlbx.Date, 10));