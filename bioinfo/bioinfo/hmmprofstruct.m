function model=hmmprofstruct(model,varargin)
%HMMPROFSTRUCT creates a profile HMM structure.
%
%   HMMPROFSTRUCT(LENGTH) returns a structure with the fields containing
%   the required parameters of a profile HMM. LENGTH specifies the number
%   of MATCH states in the model. All other mandatory model parameters are
%   initialized to default values.
%
%   HMMPROFSTRUCT(LENGTH,'field1',VALUE1,'field2',VALUE2,...) creates a
%   profile HMM using the specified fields and parameters. All other
%   mandatory model parameters are initialized to default values.
%
%   HMMPROFSTRUCT(MODEL,'field1',VALUES1,'field2',VALUES2,...) returns the
%   updated profile HMM with the specified fields and parameters updated.
%   All other mandatory model parameters are taken from the reference
%   MODEL. 
%
%   HMM PROFILE STRUCTURE FORMAT:
%
%   1. Model parameters fields (mandatory)
%      All probability values are in the range [0 1]. Generally
%      distributions should add up to 1, if not HMMPROFSTRUCT adjusts them.
%
%   ModelLength: Length of the profile (number of MATCH states).
%   Alphabet: 'AA' or 'NT',  default='AA'.
%   MatchEmission: Symbol emission probabilities for the MATCH states, size
%        is [ModelLength x AlphaLength]. Defaults to uniform distributions.
%        May accept a structure with residue counts (see AACOUNT or
%        BASECOUNT).  
%   InsertEmission: Symbol emission probabilities in the INSERT states,
%        size is [ModelLength x AlphaLength]. Defaults to uniform
%        distributions. May accept a structure with residue counts (see 
%        AACOUNT or BASECOUNT). 
%   NullEmission: Symbol emission probabilities in the MATCH and INSERT
%        states for the NULL model, size is [1 x AlphaLength]. Defaults to
%        a uniform distribution. May accept a structure with residue counts
%        (see AACOUNT or BASECOUNT). The NULL model is used to compute the
%        log-odds ratio at every state and avoid overflow when propagating
%        the probabilities through the model. 
%   BeginX: BEGIN state transition probabilities, format is
%        [B->D1 B->M1 B->M2 B->M3 .... B->Mend]. Defaults to 
%        [0.01 0.99 0 0 ... 0]. 
%   MatchX: MATCH state transition probabilities, format is:
%                  [ M1->M2 M2->M3 ... M[end-1]->Mend      ;
%                    M1->I1 M2->I2 ... M[end-1]->I[end-1]  ;
%                    M1->D2 M2->D3 ... M[end-1]->Dend      ;
%                    M1->E  M2->E  ... M[end-1]->E         ]
%        Defaults to repmat([0.998 0.001 0.001 0],profLength-1,1).
%   InsertX: INSERT state transition probabilities, format is:
%                   [ I1->M2 I2->M3 ... I[end-1]->Mend     ;
%                     I1->I1 I2->I2 ... I[end-1]->I[end-1] ]
%        Defaults to repmat([0.5 0.5],profLength-1,1).
%   DeleteX: DELETE state transition probabilities format is: 
%                   [ D1->M2 D2->M3 ... D[end-1]->Mend  ;
%                     D1->D2 D2->D3 ... D[end-1]->Dend  ]
%        Defaults to repmat([0.5 0.5],profLength-1,1).
%   FlankingInsertX: Flanking insert states (N and C) used for LOCAL
%        profile alignment, format is [N->B C->T; N->N C->C]. Defaults to 
%        [0.01 0.01; 0.99 0.99], to force global alignment use 
%        S.FlankingInsertsX = [1 1; 0 0].
%   LoopX: Loop states transition probabilities used for multiple hits
%        alignment, format is [E->C J->B; E->J J->J]. Defaults to 
%        [0.5 0.01; 0.5 0.99].
%   NullX: Null transition probabilities used to provide scores with
%        log-odds values also for state transitions, Format is 
%        [G->F ; G->G]. Defaults to [0.01; 0.99].
%
%   2. Annotation fields (optional)
%
%   Name: Model name.
%   IDNumber: Identification Number.
%   Description: Short description of the model.
%
%   Example:
%      hmmprofstruct(100,'Alphabet','AA')
%
%   See also GETHMMPROF, HMMPROFALIGN, HMMPROFESTIMATE, HMMPROFGENERATE,
%   HMMPROFMERGE, PFAMHMMREAD, SHOWHMMPROF, AACOUNT, BASECOUNT.

% Copyright 2003-2008 The MathWorks, Inc.


if (nargin==0) && (nargout==0)
    help hmmprofstruct %#ok<MCHLP>
    return
end

% defining some constants
fieldKeys = {'ModelLength','Alphabet','MatchEmission','InsertEmission','NullEmission',...
    'BeginX','MatchX','InsertX','DeleteX','FlankingInsertX','LoopX','NullX',...
    'XB','XM','XI','XD','XF','XL','XN'};

% check varargin
if rem(nargin,2) == 0
    error(message('bioinfo:hmmprofstruct:IncorrectNumberOfArguments', mfilename));
end
if nargin>1 && ~ischar([varargin{1:2:end}])
    error(message('bioinfo:hmmprofstruct:IncorrectParameterNames', mfilename));
end

% first input, reference model or model length ?
if isstruct(model)
    % Validate HMM model:
    model = checkhmmprof(model);

    profLength = model.ModelLength;
    referenceModel = true;
elseif  isnumeric(model) && numel(model)==1 && ismember(model,2:10000)
    profLength = model;
    referenceModel = false;
    model = struct; % initialize to empty struct
else
    error(message('bioinfo:hmmprofstruct:NoValidInput'));
end

% Now setting the alphabet
if referenceModel % take from ref and ignore any 'alphabet' varargin
    isamino = strcmpi(model.Alphabet,'AA');
else
    isamino = true; %default value 'alphabet' varargin
    for j=1:2:nargin-1 % look for
        pname = varargin{j};
        if strmatch(lower(pname),lower(fieldKeys))==2  
            isamino = strcmpi(varargin{j+1},'AA');
        end
    end
end
alphaLength = 4 + isamino * 16;

% Set the default model (if no reference was given)
if ~referenceModel
    model.ModelLength = profLength;
    if isamino
        model.Alphabet = 'AA';
    else
        model.Alphabet = 'NT';
    end
    model.MatchEmission = ones(profLength,alphaLength)/alphaLength;
    model.InsertEmission = ones(profLength,alphaLength)/alphaLength;
    model.NullEmission = ones(1,alphaLength)/alphaLength;
    model.BeginX = [0.01;0.99;zeros(profLength-1,1)];
    model.MatchX = repmat([0.998 0.001 0.001 0],profLength-1,1);
    model.InsertX = repmat([0.5 0.5],profLength-1,1);
    model.DeleteX = repmat([0.5 0.5],profLength-1,1);
    model.FlankingInsertX = [0.01 0.01; 0.99 0.99];
    model.LoopX = [0.5 0.01; 0.5 0.99];
    model.NullX = [0.01; 0.99];
end

% now updating structure depending on the remaining varargin
for j=1:2:nargin-2
    pname = varargin{j};
    pval = varargin{j+1};
    k = find(strncmpi(pname, fieldKeys,numel(pname))); 
    if isempty(k)
        error(message('bioinfo:hmmprofstruct:UnknownParameterName', pname));
    elseif length(k)>1
        error(message('bioinfo:hmmprofstruct:AmbiguousParameterName', pname));
    else % length(k)==1
        switch(k)
            case 1 %'ModelLength'
                error(message('bioinfo:hmmprofstruct:ModelLengthNotAllowed'));
            case 2 %'Alphabet'
                if strcmpi(pval,'AA')~=isamino
                    error(message('bioinfo:hmmprofstruct:ChangeAlphaNotAllowed'));
                    %else ignore, already taken care of
                end
            otherwise
                % convenient to set the distributions using aacount or basecount
                if isstruct(pval)
                    pval=cell2mat(struct2cell(pval));
                end
                if ~isnumeric(pval)
                    error(message('bioinfo:hmmprofstruct:NonNumericParameter', fieldKeys{ k }))
                end
                if ~all(pval(:)>=0 & pval(:)<inf)
                    error(message('bioinfo:hmmprofstruct:OutOfRangeParameter', fieldKeys{ k }))
                end
                switch(k)
                    case 3 %'MatchEmission'
                        if numel(pval)== alphaLength
                            model.MatchEmission=repmat(pval(:)'/sum(pval),profLength,1);
                        elseif all(size(pval)==[profLength alphaLength])
                            pval=diag(1./sum(pval,2))*pval;
                            model.MatchEmission=pval;
                        else
                            error(message('bioinfo:hmmprofstruct:InconsistentSizeMatchEmission', fieldKeys{ k }))
                        end
                    case 4 %'InsertEmission'
                        if numel(pval)== alphaLength
                            model.InsertEmission=repmat(pval(:)'/sum(pval),profLength,1);
                        elseif all(size(pval)==[profLength alphaLength])
                            pval=diag(1./sum(pval,2))*pval;
                            model.InsertEmission=pval;
                        else
                            error(message('bioinfo:hmmprofstruct:InconsistentSizeInsertEmission', fieldKeys{ k }))
                        end
                    case 5 %'NullEmission'
                        if numel(pval)== alphaLength
                            model.NullEmission=pval(:)'/sum(pval);
                        else
                            error(message('bioinfo:hmmprofstruct:InconsistentSizeNullEmission', fieldKeys{ k }))
                        end
                    case {6,13} %'BeginX'
                        if numel(pval)== (profLength+1)
                            model.BeginX=pval(:)/sum(pval);
                        elseif numel(pval)==2
                            model.BeginX=[pval(:)/sum(pval);zeros(profLength-1,1)];
                        else
                            error(message('bioinfo:hmmprofstruct:InconsistentSizeBeginX', fieldKeys{ k }))
                        end
                    case {7,14} %'MatchX'
                        if numel(pval)==4
                            model.MatchX=repmat(pval(:)'/sum(pval),profLength-1,1);
                        elseif numel(pval)==3
                            model.MatchX=repmat([pval(:)' 0]/sum(pval),profLength-1,1);
                        elseif all(size(pval)==[profLength-1 4])
                            pval=diag(1./sum(pval,2))*pval;
                            model.MatchX=pval;
                        elseif all(size(pval)==[profLength-1 3])
                            pval=diag(1./sum(pval,2))*pval;
                            model.MatchX=[pval;zeros(profLength-1,1)];
                        else
                            error(message('bioinfo:hmmprofstruct:InconsistentSizeMatchX', fieldKeys{ k }))
                        end
                    case {8,15} %'InsertX'
                        if numel(pval)==2
                            model.InsertX=repmat(pval(:)'/sum(pval),profLength-1,1);
                        elseif all(size(pval)==[profLength-1 2])
                            pval=diag(1./sum(pval,2))*pval;
                            model.InsertX=pval;
                        else
                            error(message('bioinfo:hmmprofstruct:InconsistentSizeInsertX', fieldKeys{ k }))
                        end
                    case {9,16} %'DeleteX'
                        if numel(pval)==2
                            model.DeleteX=repmat(pval(:)'/sum(pval),profLength-1,1);
                        elseif all(size(pval)==[profLength-1 2])
                            pval=diag(1./sum(pval,2))*pval;
                            model.DeleteX=pval;
                        else
                            error(message('bioinfo:hmmprofstruct:InconsistentSizeDeleteX', fieldKeys{ k }))
                        end
                    case {10,17} %'FlankingInsertX'
                        if numel(pval)==2
                            model.FlankingInsertX=repmat(pval(:)/sum(pval),1,2);
                        elseif all(size(pval)==[2 2])
                            pval=pval*diag(1./sum(pval));
                            model.FlankingInsertX=pval;
                        else
                            error(message('bioinfo:hmmprofstruct:InconsistentSizeFlankingInsertX', fieldKeys{ k }))
                        end
                    case {11,18} %'LoopX'
                        if numel(pval)==2
                            model.LoopX=repmat(pval(:)/sum(pval),1,2);
                        elseif all(size(pval)==[2 2])
                            pval=pval*diag(1./sum(pval));
                            model.LoopX=pval;
                        else
                            error(message('bioinfo:hmmprofstruct:InconsistentSizeLoopX', fieldKeys{ k }))
                        end
                    case {12,19} %'NullX'
                        if numel(pval)==2
                            model.NullX=pval(:)/sum(pval);
                        else
                            error(message('bioinfo:hmmprofstruct:InconsistentSizeNullX', fieldKeys{ k }))
                        end
                end %second switch
        end % first switch
    end %if length(k)==1
end % for all varargin
