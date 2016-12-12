function openvar(name, tr) %#ok
%OPENVAR Opens a phylogenetic tree object for graphical editing.

% Copyright 2003-2006 The MathWorks, Inc.


try
    view(tr);
catch theException
    % rethrows the error into a dlg window
    errordlg(theException.message, 'Inspection error', 'modal');
end
