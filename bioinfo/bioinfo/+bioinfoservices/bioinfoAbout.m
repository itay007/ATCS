function aboutstr = bioinfoAbout(varargin)%#ok
%BIOINFOABOUT Return a string about the Bioinformatics Toolbox.

tlbx = ver('bioinfo');
tlbx = tlbx(1);
aboutstr = sprintf('%s\nVersion %s\nCopyright 2003-%s, The MathWorks, Inc.', ...
    tlbx.Name, tlbx.Version, datestr(tlbx.Date, 10));