function model = checkhmmprof(model)
%CHECKHMMPROF Validates a HMM profile model.
%   
%   CHECKHMMPROF(STRUCTURE) returns a valid HMM profile model or errors out 
%   if the model cannot be validated.
%
%   The function actualizes a HMM profile model saved on R13 into R14 valid 
%   models by initializing the fields LoopX and NullX to their default
%   values if they are not present in the input structure.

%   Copyright 2003-2012 The MathWorks, Inc.


% model must be an structure
if ~isstruct(model)
    msgId = 'bioinfo:checkhmmprof:ModelIsNotAStructure';
    x = MException(msgId,getString(message(msgId)));
    x.throwAsCaller;
end

% not allowing arrays of structures
if numel(model)>1
    model = model(1);
    warning(message('bioinfo:checkhmmprof:NoMultipleModels',upper(evalin('caller','mfilename'))))
end    

% look for the minimum required fields
fieldKeys = {'ModelLength','Alphabet','MatchEmission','InsertEmission',...
             'NullEmission','BeginX','MatchX','InsertX','DeleteX',...
             'FlankingInsertX','LoopX','NullX'};
missingFields = ~ismember(fieldKeys,fieldnames(model));
if any(missingFields(1:10))
    error(message('bioinfo:checkhmmprof:NoValidModel'))
end

% checking missing fields not used in R13 and initializing to defaults
if missingFields(11)
    warning(message('bioinfo:checkhmmprof:LoopXMissingInModel'))
    model.LoopX = [0.5 0.01; 0.5 0.99];
end
if missingFields(12)
    warning(message('bioinfo:checkhmmprof:NullXMissingInModel'))
    model.NullX = [0.01; 0.99];
end
