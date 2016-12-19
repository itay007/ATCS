function scaleimagefigure(hfig, hax, himg)
%SCALEIMAGEFIGURE resizes the figure window to fit the image size.
%
%	SCALEIMAGEFIGURE resize the figure HFIG and the axis HAX to fit the image
%	HIMG dimension.
%   
%   SCALEIMAGEFIGURE is a helper function for CONNECTKEGG example.

%   Copyright 2006-2012 The MathWorks, Inc. 


%  Set an edge around the image. 
edge = 2;
set(hax, 'DataAspectRatio', [1 1 1]);
%%
% Get the image dimension.
img_width  = size(get(himg, 'CData'),2);
img_height = size(get(himg, 'CData'),1);

% Set the figure dimension
fig_width = img_width + edge;
fig_height = img_height + edge;

% Get the current figure and axis positions
fig_unit = get(hfig, 'Units');
ax_unit = get(hax, 'Units');

set(hfig, 'Units', 'pixels')
set(hax, 'Units', 'pixels')
fpos = get(hfig, 'Position');
apos = get(hax, 'position');

% Determine the new figure and axis positions
new_fpos = fpos;
new_fpos(1) = fpos(1) - (fig_width - fpos(3))/2;
new_fpos(2) = fpos(2) - (fig_height - fpos(4))/2;
new_fpos(3) = fig_width;
new_fpos(4) = fig_height;

apos = [edge/2, edge/2, img_width, img_height];
%%
% Depends on you screen size, you can translate the figure position if
% necessary. Get the screen dimension
win_unit = get(0, 'Units');
set(0, 'Units', 'pixels');
screen_size = get(0, 'ScreenSize');
screen_width = screen_size(3);
screen_height = screen_size(4);

% Check the width
side = max(edge, 100);
delta = (screen_width - side) - (new_fpos(1) + new_fpos(3));
if (delta < 0)
    new_fpos(1) = new_fpos(1) + delta;
end

% Check the height
delta = (screen_height-side) - (new_fpos(2) + new_fpos(4));
if (delta < 0)
    new_fpos(2) = new_fpos(2) + delta;
end

%%
% Set the figure and axis positions
set(hfig, 'Position', new_fpos)
set(hax, 'Position', apos)

set(hfig, 'Units', fig_unit);
set(hax, 'Units', ax_unit);
set(0, 'Units', win_unit);
