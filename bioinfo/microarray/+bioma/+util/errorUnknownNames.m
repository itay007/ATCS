function errorUnknownNames(ukNames)
%ERRORUNKNOWNNAMES Throw error message with listed names.
%
% ERRORUNKNOWNNAMES(NAMES) Throw error of message list the unrecognized
% NAMES.

%   Copyright 2009 The MathWorks, Inc. 


ukn = ''; % list with missing fields
if ischar(ukNames)
   ukNames = {ukNames}; 
end
    
for i=1:numel(ukNames)
    ukn = [ukn ukNames{i} ', ']; %#ok<AGROW>
end
ukn = ukn(1:end-2); % remove extra comma and space from the end
    error(message('bioinfo:errorUnknownNames:UnrecognizableNames', ukn));
end
