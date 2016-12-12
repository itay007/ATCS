function disableHGBehaviors(h)
%DISABLEHGBEHAVIORS disables hg object behaviors: datacursor, legend

% Copyright 2009 The MathWorks, Inc.


%== Turn off data cursor
hb=hggetbehavior(h, 'Datacursor');
set(hb, 'Enable', false)

%== Turn off legend
set(get(get(h,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off');
%== Turn data brushing
end
