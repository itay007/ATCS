function fixGenomicPositionLabels(ha)
% FIXGENOMICPOSITIONLABELS Helper function that adds callback functions for
% updating XTickLabel and datacursor, so the X axis labels are formatted
% properly as genomic positions. FIXGENOMICPOSITIONLABELS also turns the
% box on and resizes the figure, so the rendering is more appealing to
% genomic data.
%
% FIXGENOMICPOSITIONLABELS(HA) HA is the axes containing the genomic data.

%   Copyright 2012 The MathWorks, Inc.

if nargin == 0
    ha = gca;
end
hf = get(ha,'Parent');   % handle to figure
dc = datacursormode(hf); % handle to datacursor manager
ad.UpdateFunctionToHonor = get(dc,'UpdateFcn'); 
setappdata(hf,'fixGenomicPositionLabels',ad)  % store old datacursor function handle
set(dc,'UpdateFcn',@updateDataCursor) 
datacursormode(hf,'on')
setappdata(ha,'datacursorHandledByFixGenomicPositionLabels',true);
lh = addlistener(ha,'XTick','PostSet',@updateXTickLabels);
setappdata(ha,'XtickListener',lh);
box(ha,'on')             % box on
set(hf,'Position',max(get(hf,'Position'),[0 0 900 0])) %resize window

end

function updateXTickLabels(~,e)
% Callback function to reformat the XTick labels

ha = get(e,'AffectedObject');
v = get(ha,'XTick');
set(ha,'XTickLabel',num2str(v'));
end

function txt = updateDataCursor(h,e)
% Callback function to prepare the text that goes in a datacursor

% find the first parent axes to the target (calling target may be inside
% hggropus)
ha = get(e,'Target');
while ~strcmp(get(ha,'Type'),'axes') && ha~=0
    ha = get(ha,'Parent');
end
% Check if fixGenomicPositionLabels owns the data cursor for this axes
if ~isempty(getappdata(ha,'datacursorHandledByFixGenomicPositionLabels'))
    % fixGenomicPositionLabels prepares the datacursor for this axes
    txt = {['Position: ',num2str(get(e,'Position')*[1;0])],...
           ['Y: ',num2str(get(e,'Position')*[0;1])]};
else
    hf = get(ha,'Parent'); %handle to figure
    ad = getappdata(hf,'fixGenomicPositionLabels');
    % check is there is another updateDataCursor that neeeds to be called
    if ~isempty(ad.UpdateFunctionToHonor)
        % Honor the previously saved UpdateFcn in axes that I do not own
        txt = feval(ad.UpdateFunctionToHonor,h,e);
    else
        % fixGenomicPositionLabels does not own this axes and the
        % previously saved UpdateFcn is empty, then just replicate MATLAB's
        % default datacursor
        pos = get(e, 'position');
        txt = {['X=',num2str(pos(1))], ['Y=',num2str(pos(2))]};
    end
end
end