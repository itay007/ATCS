function parentResizeCB(src, evt, hcom, varargin) %#ok
%PARENTRESIZECB Callback when parent of an UI componet is resized.
% 
%   PARENTRESIZECB(SRC, EVT, HCOM, WIDTHRATIO, HEIGHTRATIO) resize UI
%   component HCOM when its parent is resized. The margin between the
%   component HCOM and the parent is set by the ratios of their width and
%   height, WIDTHRATIO, HEIGHTRATIO. The default ratio is 1/25. 

%   Copyright 2009 The MathWorks, Inc. 


hpnt = get(hcom, 'Parent');
%== Advoid warning incase this property was not set to 'Auto'
if strcmpi(get(hpnt, 'Type'), 'figure')
    set(hpnt, 'PaperPositionMode', 'auto');
end

pPos = get(hpnt, 'Position');
widthRatio = 1/25;
heightRatio = 1/25;
if nargin >= 4
    widthRatio = varargin{1};
end

if nargin == 5
    heightRatio = varargin{2};
end

set(hcom, 'position', [pPos(3)*widthRatio
                       pPos(4)*heightRatio 
                       pPos(3)*(1-2*widthRatio)
                       pPos(4)*(1-2*heightRatio)]);

end
