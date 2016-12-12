function schema
%

% Copyright 2003-2004 The MathWorks, Inc.

schema.package('biograph');

%% set shape types
if (isempty(findtype('BiographShape')))
    schema.EnumType('BiographShape',...
        {'box','ellipse','circle','rect','rectangle','diamond',...
         'trapezium','house','invtrapezium','invhouse','parallelogram'});
% Future shapes:     
%          'polygon','egg','triangle','invtriangle'
%          'pentagon','hexagon','septagon','octagon'
end

