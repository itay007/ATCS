function hg_handles = customnodedraw(node_handle)
%CUSTOMNODEDRAW User function example to draw customized nodes.
%
%  Utilize a user specified function to customize the node drawing of a 
%  BIOGRAPH object. The function must have the following form:
%
%     function  HG_HANDLES = CUSTOMNODEDRAW(NODE_HANDLE)
%
%  taking as input argument a @biograph/@node object. The function
%  CUSTOMNODEDRAW uses the @biograph/@node properties for drawing the node,
%  such as the "Position", the "Size" and the "Label" of the node. The user
%  can store application data into the property "UserData" before creating
%  the graph layout, then CUSTOMNODEDRAW can utilize it to enhance the node
%  rendered information, for example confidence values, rates or
%  probability distributions.
%
%  The handles of all objects drawn by CUSTOMNODEDRAW must be returned into
%  the vector HG_HANDLES, so the @biograph object deletes them when
%  necessary.


%   Copyright 2003-2008 The MathWorks, Inc.



% Get the axes handle where all hg objects are attached, all hg objects
% should specifically set "haxes" as the parent to ensure that the
% customized node is rendered in the proper figure.
haxes = node_handle.up.hgAxes;

% Get the current scale between the computed layout and rendered graph
biographScale = node_handle.up.Scale;

% Get the size of the node
rx = node_handle.Size(1);
ry = node_handle.Size(2);

% Get the center of the node
x = node_handle.Position(1)*biographScale;
y = node_handle.Position(2)*biographScale;

% In this example you create a pie chart with the vector stored in
% node_handle.UserData.Distribution

if isempty(node_handle.UserData) || ...
        ~isfield(node_handle.UserData,'Distribution') || ...
        isempty(node_handle.UserData.Distribution)

    % could not find data to draw the pie, then just draw a circle
    t = [0:pi/25:2*pi];
    px = [rx/2*cos(t)+x];
    py = [ry/2*sin(t)+y];
    hg_handles = zeros(2,1); % initialize output handles
    hg_handles(1) = patch(px,py,[1 0 0],'Parent',haxes);
    
else % draw the pie chart
    
    d = node_handle.UserData.Distribution;
    d = d(:)./sum(d(:)); % normalize

    alpha = [0;cumsum(d)]*2*pi;

    % Draw each slice of the pie chart
    num_slices = numel(d);
    colors = jet(num_slices);
    hg_handles = zeros(num_slices+1,1); % initialize output handles

    for i = 1:num_slices
        t = [alpha(i):pi/25:alpha(i+1) alpha(i+1)];
        px = [x rx/2*cos(t)+x x];
        py = [y ry/2*sin(t)+y y];
        hg_handles(i) = patch(px,py,colors(i,:),'Parent',haxes);
    end

end


% Place the label in the lower-right corner of the node
hg_handles(end) = text(x+rx/2,y-ry/2,node_handle.id,...
    'HorizontalAlignment','Left',...
    'VerticalAlignment','Middle',...
    'Fontsize',7,'Interpreter','none',...
    'Parent',haxes);


