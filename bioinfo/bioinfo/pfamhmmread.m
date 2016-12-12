function data=pfamhmmread(filename)
%PFAMHMMREAD reads in PFAM-HHM format data files.
%
%   S = PFAMHMMREAD(FILE) reads in a Hidden Markov Model described by the
%   PFAM format from FILE, and converts it to structure S, containing
%   fields corresponding to annotations and parameters of the model. For
%   more information about the model structure format see HMMPROFSTRUCT.
%   FILE can also be an URL or a MATLAB cell array that contains the text
%   of a PFAM format file.
%
%   If a file contains multiple models, they are loaded into an array of
%   structures. However you should not try to load files with more than
%   1000 models.
%
%   Examples:
%
%       pfamhmmread('pf00002.ls')
%       site='http://pfam.sanger.ac.uk/';
%       pfamhmmread([site '/family/hmm?mode=ls&id=7tm_2'])
%
%   See also GETHMMALIGNMENT, GETHMMPROF, HMMPROFALIGN, HMMPROFDEMO,
%   HMMPROFSTRUCT, SHOWHMMPROF.

%   PFAMHMMREAD is based on the HMMER2.0 and HMMER3/b file formats.

% Copyright 2003-2012 The MathWorks, Inc.


% check input is char
if ~ischar(filename)
    error(message('bioinfo:pfamhmmread:InvalidInput'))
end

% See if we have a URL
if size(filename,1)>1  % is padded string
    pftext = cell(size(filename,1),1);
    for i=1:size(filename,1)
        pftext(i)=strtrim(strread(filename(i,:),'%s','whitespace',''));
    end
    % try then if it is an url
elseif (strfind(filename(1:min(10,end)), '://'))
    if (~usejava('jvm'))
        error(message('bioinfo:pfamhmmread:NoJava'))
    end
    try
        pftext = urlread(filename);
    catch theException
        error(message('bioinfo:pfamhmmread:CannotReadURL', filename));
    end
    pftext = textscan(pftext,'%s','delimiter','\n');
    pftext = pftext{1};
    % try then if it is a valid filename
elseif  (exist(filename,'file') || exist(fullfile(pwd,filename),'file'))
    fid = fopen(filename,'r');
    pftext = textscan(fid,'%s','delimiter','\n');
    pftext = pftext{1};
    fclose(fid);
else % must be a string with '\n', convert to cell
    pftext = textscan(filename,'%s','delimiter','\n');
    pftext = pftext{1};
end

% At this point pftext must be a cellstr
if ~iscellstr(pftext) 
    error(message('bioinfo:pfamhmmread:PfamNotValid'))
end

% finding closing markers of every record
recordLimits =[0;strmatch('//',pftext,'exact')];
numberRecords = length(recordLimits)-1;

if numberRecords==0
    error(message('bioinfo:pfamhmmread:NoRecordsFound'))
end

try
    data(numberRecords) = struct;
    for count = 1:numberRecords
        pf = pftext(recordLimits(count)+1:recordLimits(count+1)-1);
        if strncmp(pf{1},'HMMER2.0',8)
            data1 = parsehmmer20(pf);
        elseif strncmp(pf{1},'HMMER3/b',8)
            data1 = parsehmmer30b(pf);
        else
            error(message('bioinfo:pfamhmmread:PfamRecordNotValid', count))
        end
        fnames = fieldnames(data1);
        for j = 1:numel(fnames)
            data(count).(fnames{j}) = data1.(fnames{j});
        end
    end %for count

catch le
    if strncmpi(le.identifier,'bioinfo:',8) % did we catch one of our errors ?
        rethrow(le)
    else
        error(message('bioinfo:pfamhmmread:IncorrectDataFormat'))
    end
end
end

function data = parsehmmer20(pf)
% Reads HMMER2.0 formatted files as the specification in "HMMER's User
% Guide", Version 2.3.2, October 2003.

        % constants in the HMMER2.0 file format
        numX = 9 ;    %number of possible transitions within a state
        scale = 1000; %scale used to recover probabilities from log-prob

        %%%%% READ ANNOTATION FIELDS %%%%%

        % Model Name -- Mandatory -- Saved in structure
        h=strmatch('NAME',pf);
        data.Name=pf{h}(7:end);

        % PFAM Accession Number -- Optional -- Saved in structure
        h=strmatch('ACC',pf);
        if isempty(h)
            data.PfamAccessionNumber=[]; 
        else
            data.PfamAccessionNumber=pf{h(1)}(7:end);
        end

        % HMM Model description -- Optional -- Saved in structure
        h=strmatch('DESC',pf);
        if isempty(h)
            data.ModelDescription=[];
        else
            data.ModelDescription=pf{h(1)}(7:end);
        end

        %%%%% READ LENGTH AND ALPHABET %%%%%

        % HMM Model length -- Mandatory -- Saved in structure
        h=strmatch('LENG',pf);
        profLength=str2num(pf{h(1)}(7:end)); 
        data.ModelLength=profLength;

        % Alphabet -- Mandatory -- Saved in structure
        data.Alphabet='AA';
        alphaLength=20;
        isAmino=true;
        h=strmatch('ALPH',pf);
        if strcmp('Nucleic',pf{h(1)}(7:end))
            data.Alphabet='NT';
            alphaLength=4;
            isAmino=false;
        end

        %%%%% READ OTHER ANNOTATIONS NOT STORED IN THE STRUCTURE %%%%%

        % Reference Annotation -- Optional -- Not used
        % h = strmatch('RF',pf);

        % Consensus Structure Annotation -- Optional -- Not used
        % h=strmatch('CS',pf);

        % Map Annotation flag -- Optional -- Used to remap but not stored
        %                                    in structure
        remap = false;
        h=strmatch('MAP',pf);
        if ~isempty(h) && strcmp('yes',pf{h(1)}(7:end))
            remap=true;
        end

        % Gathering cutoff -- Optional -- Not used
        h=strmatch('GA',pf);
        if isempty(h)
            % GatheringCutoff=zeros(0,2); % Not used
        else
            GatheringCutoff=str2num(pf{h(1)}(7:end)); %#ok<*ST2NM>
            if ~all(size(GatheringCutoff)==[1 2])
                error(message('bioinfo:pfamhmmread:IncorrectDataFormatGC'))
            end
        end

        % Trusted cutoff -- Optional -- Not used
        h=strmatch('TC',pf);
        if isempty(h)
            % TrustedCutoff=zeros(0,2); % Not uesed
        else
            TrustedCutoff=str2num(pf{h(1)}(7:end));
            if ~all(size(TrustedCutoff)==[1 2])
                error(message('bioinfo:pfamhmmread:IncorrectDataFormatTC'))
            end
        end

        % Noise cutoff -- Optional -- Not used
        h=strmatch('TC',pf);
        if isempty(h)
            % NoiseCutoff=zeros(0,2); % Not used
        else
            NoiseCutoff=str2num(pf{h(1)}(7:end));
            if ~all(size(NoiseCutoff)==[1 2])
                error(message('bioinfo:pfamhmmread:IncorrectDataFormatNC'))
            end
        end

        %%%%% READ PROBABILITIES %%%%%

        % Special Transition Log-Probabilities -- Mandatory
        h=strmatch('XT',pf);
        XSpecial=str2num(pf{h(1)}(7:end));
        if ~all(size(XSpecial)==[1 8])
            error(message('bioinfo:pfamhmmread:IncorrectDataFormatXT'))
        end

        % Null Model Transition Log-Probabilities -- Mandatory
        h=strmatch('NULT',pf);
        NullTran=str2num(pf{h(1)}(7:end));
        if ~all(size(NullTran)==[1 2])
            error(message('bioinfo:pfamhmmread:IncorrectDataFormatNULT'))
        end

        % Null Model Symbol Emission Log-Probabilities -- Mandatory
        % reorder at the end
        h=strmatch('NULE',pf);
        NullEmission=str2num(pf{h(1)}(7:end));
        if strcmp(data.Alphabet,'NT') && ~all(size(NullEmission)==[1 4])
            error(message('bioinfo:pfamhmmread:IncorrectDataFormatNULENT'))
        end
        if strcmp(data.Alphabet,'AA') && ~all(size(NullEmission)==[1 20])
            error(message('bioinfo:pfamhmmread:IncorrectDataFormatNULEAA'))
        end

        % Extreme value distribution parameter -- Optional
        h=strmatch('EVD',pf);
        if isempty(h)
           % ExtremeValue=zeros(0,2); % Not used
        else
            ExtremeValue=str2num(pf{h(1)}(7:end));
            if ~all(size(ExtremeValue)==[1 2])
                error(message('bioinfo:pfamhmmread:IncorrectDataFormatEVD'))
            end
        end

        % Looks for the HMM flag -- Mandatory
        h=strmatch('HMM',pf);
        if length(h)~=2 % it must have two, one for HMMER2.0 and the other the flag
            error(message('bioinfo:pfamhmmread:IncorrectDataFormatHMM'))
        end
        h = h(2); % using the HMM flag;

        % Alphabet order, re-order as aa2int (or nt2int) -- Mandatory
        % reorderAlpha is used later to read in the Symbol Emission lines
        if isAmino
            reorderAlpha=aa2int(removeallblanks(pf{h}(7:end)));
        else
            reorderAlpha=nt2int(removeallblanks(pf{h}(7:end)));
        end

        % Identify the state transitions in the header
        % reorderStateX is used later to read in the State Transition lines
        StateAlpha = ('MIDBE')';
        StateAlphaI = zeros(77,1);
        StateAlphaI('MIDBE') = 1:5;
        %StateXMap =[1 1;1 2;1 3;1 5;2 1;2 2;3 1;3 3;4 1];
        temp = upper(pf{h+1}(any(repmat(upper(pf{h+1}),5,1)==repmat(StateAlpha,1,size(pf{h+1},2)))));
        [~,reorderStateX] = sortrows(StateAlphaI(reshape(temp,2,length(temp)/2)'));

        % Transitions of the B state -- Mandatory
        XBState = CorrectStars(pf{h+2});
        h=h+3;

        % pre-allocate
        MatchEmission = zeros(profLength,alphaLength);
        InsertEmission = zeros(profLength,alphaLength);
        XState = zeros(profLength,numX);
        Map = zeros(profLength,1);

        for state = 1:profLength
            % Every state of the model is described by three lines

            % Match emission line
            p = CorrectStars(pf{h});
            MatchEmission(state,reorderAlpha) = p(2:alphaLength+1);
            if remap
                Map(state) = p(alphaLength+2);
            else
                Map(state) = p(1);
            end

            % Insert emission line
            p = CorrectStars(pf{h+1}(2:end)); % after the '-' or the RF
            InsertEmission(state,reorderAlpha) = p;

            % State transition line
            p=CorrectStars(pf{h+2}(2:end)); % after the '-' or the CS
            XState(state,:)=p(reorderStateX);

            h=h+3;

        end %for state

        NullEmission(reorderAlpha)=NullEmission;

        %%%RE-SCALE ALL PROBABILITIES AND UNDO LOGS

        %%% NULL MODEL EMISSION PROB
        nep=(2.^(NullEmission/scale))/alphaLength;
        nep=nep/sum(nep);              %re-adjusting to sum one

        %%% MATCH EMISSION PROB
        mep=(2.^(MatchEmission/scale))*diag(nep);
        mep=diag(1./sum(mep,2))*mep;   %re-adjusting to sum one

        %%% INSERT EMISSION PROB
        iep=(2.^(InsertEmission/scale))*diag(nep);
        iep(~sum(iep,2),:)=1/alphaLength; %correct rows that sum zero
        iep=diag(1./sum(iep,2))*iep;   %re-adjusting to sum one

        %%% NULL MODEL TRANSITION PROB
        nxp=2.^(NullTran([2,1])/scale);
        nxp=nxp/sum(nxp);             %re-adjusting to sum one
        
        %%% BEGIN-[DEL_1-MATCH_X] TRANSITION PROB
        %%% [B->D1, B->M1, B->M2, B->M3, B->M4 ... ]
        bxp=[2.^(XBState(3)/scale);2.^(XState(:,9)/scale)];
        bxp=bxp/sum(bxp);             %re-adjusting to sum one

        %%% MATCH TRANSITION PROB
        mxp=2.^(XState(:,[1 2 3 4])/scale);
        mxp(~sum(mxp,2),:)=1/4;  %correct rows that sum zero
        mxp=diag(1./sum(mxp,2))*mxp;   %re-adjusting to sum one

        %%% INSERT TRANSITION PROB
        ixp=2.^(XState(:,[5 6])/scale);
        ixp(~sum(ixp,2),:)=1/2;  %correct rows that sum zero
        ixp=diag(1./sum(ixp,2))*ixp;   %re-adjusting to sum one

        %%% DELETE TRANSITION PROB
        dxp=2.^(XState(:,[7 8])/scale);
        dxp(~sum(dxp,2),:)=1/2;  %correct rows that sum zero
        dxp=diag(1./sum(dxp,2))*dxp;   %re-adjusting to sum one

        %%% LEFT FLANKING INSERT TRANSITION PROB
        lixp=2.^(XSpecial([1,2])/scale);
        lixp=lixp/sum(lixp);          %re-adjusting to sum one

        %%% RIGTH FLANKING INSERT TRANSITION PROB
        rixp=2.^(XSpecial([5,6])/scale);
        rixp=lixp/sum(rixp);          %re-adjusting to sum one
        
        %%% LOOP INSERT TRANSITION PROB
        lxp=2.^(XSpecial([7,8])/scale);
        lxp=lxp/sum(lxp);          %re-adjusting to sum one
        
        %%% END TRANSITION PROB
        %%% [E->C, E->J]
        endxp=2.^(XSpecial([3,4])/scale);
        endxp=endxp/sum(endxp);          %re-adjusting to sum one

        %%% Store data in structure
        data.MatchEmission = mep;
        data.InsertEmission = iep;
        data.NullEmission = nep;
        data.BeginX = bxp;
        data.MatchX = mxp(1:end-1,:); %last state is always [0 0 0 1]
        data.InsertX = ixp(1:end-1,:); %last state does not exist
        data.DeleteX = dxp(1:end-1,:); %last state always goes to E
        data.FlankingInsertX = [lixp;rixp]';
        data.LoopX = [endxp;lxp]';
        data.NullX = nxp';
end

function data = parsehmmer30b(pf)
% Reads HMMER3/b formatted files as the specification in "HMMER3 beta test:
% User Guide", Version 3.0b3, November 2009.

        %%%%% READ ANNOTATION FIELDS %%%%%

        % Model Name -- Mandatory -- Saved in structure
        h=strmatch('NAME',pf);
        data.Name=pf{h}(7:end);

        % PFAM Accession Number -- Optional -- Saved in structure
        h=strmatch('ACC',pf);
        if isempty(h)
            data.PfamAccessionNumber=[]; 
        else
            data.PfamAccessionNumber=pf{h(1)}(7:end);
        end

        % HMM Model description -- Optional -- Saved in structure
        h=strmatch('DESC',pf);
        if isempty(h)
            data.ModelDescription=[];
        else
            data.ModelDescription=pf{h(1)}(7:end);
        end

        %%%%% READ LENGTH AND ALPHABET %%%%%

        % HMM Model length -- Mandatory -- Saved in structure
        h=strmatch('LENG',pf);
        profLength=str2num(pf{h(1)}(7:end)); 
        data.ModelLength=profLength;

        % Alphabet -- Mandatory -- Saved in structure
        data.Alphabet='AA';
        alphaLength=20;
        isAmino=true;
        h=strmatch('ALPH',pf);
        if strcmp('Nucleic',pf{h(1)}(7:end))
            data.Alphabet='NT';
            alphaLength=4;
            isAmino=false;
        end

        %%%%% READ OTHER ANNOTATIONS NOT STORED IN THE STRUCTURE %%%%%

        % Reference Annotation -- Optional -- Not used
        % h = strmatch('RF',pf);

        % Consensus Structure Annotation -- Optional -- Not used
        % h=strmatch('CS',pf);

        % Map Annotation flag -- Optional -- Used to remap but not stored
        %                                    in structure
        remap = false;
        h=strmatch('MAP',pf);
        if ~isempty(h) && strcmp('yes',pf{h(1)}(7:end))
            remap=true;
        end
        
        % Date the model was constructed (DATE) -- Optional -- Not used
        
        % Command line log (COM or BM) -- Optional -- Not used
        
        % Sequence number (NSEQ) -- Optional -- Not used
        
        % Effective sequence number (EFFN) -- Optional -- Not used

        % Gathering cutoff -- Optional -- Not used
        h=strmatch('GA',pf);
        if isempty(h)
            % GatheringCutoff=zeros(0,2); % Not used
        else
            GatheringCutoff=str2num(pf{h(1)}(7:end)); %#ok<*ST2NM>
            if ~all(size(GatheringCutoff)==[1 2])
                error(message('bioinfo:pfamhmmread:IncorrectDataFormatGC'))
            end
        end

        % Trusted cutoff -- Optional -- Not used
        h=strmatch('TC',pf);
        if isempty(h)
            % TrustedCutoff=zeros(0,2); % Not uesed
        else
            TrustedCutoff=str2num(pf{h(1)}(7:end));
            if ~all(size(TrustedCutoff)==[1 2])
                error(message('bioinfo:pfamhmmread:IncorrectDataFormatTC'))
            end
        end

        % Noise cutoff -- Optional -- Not used
        h=strmatch('TC',pf);
        if isempty(h)
            % NoiseCutoff=zeros(0,2); % Not used
        else
            NoiseCutoff=str2num(pf{h(1)}(7:end));
            if ~all(size(NoiseCutoff)==[1 2])
                error(message('bioinfo:pfamhmmread:IncorrectDataFormatNC'))
            end
        end

        % Statistical parameters needed for E-value calculations (STATS) -- Optional -- Not used      
        
        %%%%% READ PROBABILITIES %%%%%
       
        % Looks for the HMM flag -- Mandatory
        h=strmatch('HMM',pf);
        if length(h)~=2 % it must have two, one for HMMERxxx and the other the flag
            error(message('bioinfo:pfamhmmread:IncorrectDataFormatHMM'))
        end
        h = h(2); % using the HMM flag;

        % Alphabet order, re-order as aa2int (or nt2int) -- Mandatory
        % reorderAlpha is used later to read in the Symbol Emission lines
        if isAmino
            reorderAlpha=aa2int(removeallblanks(pf{h}(7:end)));
        else
            reorderAlpha=nt2int(removeallblanks(pf{h}(7:end)));
        end

        % Identify the state transitions in the header
        % reorderStateX is used later to read in the State Transition lines
        StateAlpha = ('MIDBE')';
        StateAlphaI = zeros(77,1);
        StateAlphaI('MIDBE') = 1:5;
        %StateXMap =[1 1;1 2;1 3;1 5;2 1;2 2;3 1;3 3;4 1];
        temp = upper(pf{h+1}(any(repmat(upper(pf{h+1}),5,1)==repmat(StateAlpha,1,size(pf{h+1},2)))));
        [~,reorderStateX] = sortrows(StateAlphaI(reshape(temp,2,length(temp)/2)'));

        h=h+2;
        
        % Read the background (null) probabilities (COMPO) --optional --
        % If the file does not contain the overall average match state
        % emission probabilities then the null emissions will be set to the
        % HMMER2.0 defaults, which are taken of the aminoacid composition
        % of SWISS-PROT 34:
        if strncmp(pf{h},'COMPO',5)
             p = CorrectStars3b(pf{h}(6:end));
             NullEmission(1,reorderAlpha) = p(1:alphaLength);
             h = h+1;
        else
             NullEmission = [2.541097704695023   2.916783476558513...   
                  3.183645141074092   2.927873831447472   4.189401700066573...
                  3.230086002171608   2.705373586487729   2.666557344376372...
                  3.775592833272285   2.830140078988520   2.339391875152078...   
                  2.822515460002360   3.739549179883168   3.226620266268808...
                  3.030459614170344   2.683192876709811   2.917476623739073...   
                  4.472898896915590   3.492788783603827   2.697748967501570];
        end
        
        % Read State 0 (or BEGIN node) emissions
        p = CorrectStars3b(pf{h});
        BeginEmission(1,reorderAlpha) = p(1:alphaLength);
        h = h+1;
        
        % Read State 0 (or BEGIN node) transitions
        p = CorrectStars3b(pf{h});
        XBState(1,1:numel(reorderStateX)) = p(reorderStateX);
        h = h+1;

        % pre-allocate
        MatchEmission = zeros(profLength,alphaLength);
        InsertEmission = zeros(profLength,alphaLength);
        XState = zeros(profLength,numel(reorderStateX));
        Map = zeros(profLength,1);

        for state = 1:profLength
            % Every state of the model is described by three lines

            % Match emission line
            p = CorrectStars3b(pf{h});
            MatchEmission(state,reorderAlpha) = p(2:alphaLength+1);
            if remap
                Map(state) = p(alphaLength+2);
            else
                Map(state) = p(1);
            end

            % Insert emission line
            p = CorrectStars3b(pf{h+1}); 
            InsertEmission(state,reorderAlpha) = p;

            % State transition line
            p = CorrectStars3b(pf{h+2}); 
            XState(state,:) = p(reorderStateX);

            h = h+3;

        end %for state    
   
        
        %%%RE-SCALE ALL PROBABILITIES AND UNDO LOGS

        %%% NULL MODEL EMISSION PROB
        nep = exp(-NullEmission);
        nep = nep/sum(nep);              %re-adjusting to sum one

        %%% MATCH EMISSION PROB
        mep = exp(-MatchEmission);
        mep = diag(1./sum(mep,2))*mep;   %re-adjusting to sum one

        %%% INSERT EMISSION PROB
        iep = exp(-InsertEmission);
        iep(~sum(iep,2),:) = 1/alphaLength; %correct rows that sum zero
        iep = diag(1./sum(iep,2))*iep;   %re-adjusting to sum one

        %%% NULL MODEL TRANSITION PROB
        % Null transitions are not stored in HMMER3/b, in HMMER2.0 the null
        % model uses an average model length of 350:
        nxp = [1 350]./351;
        
        %%% BEGIN-[DEL_1-MATCH_X] TRANSITION PROB
        %%% [B->D1, B->M1, B->M2, B->M3, B->M4 ... ]
        % The Bioinformatics Toolbox implements Plan-7 proposed in
        % HMMER2.0, therefore M0->M1 and M0->D1 from HMMER3/b are
        % marginalized to the Plan-7 B->M1 and B->D1 transition
        % probabilities. All B->Mx where x>1 probabilities are set to 0.
        bxp = zeros(1,profLength+1);
        bxp([1 2]) = exp(-XBState([3,1]));
        bxp = bxp/sum(bxp);             %re-adjusting to sum one

        %%% MATCH TRANSITION PROB
        % In HMMER3/b there are not Mx->E transitions
        mxp = exp(-XState(:,[1 2 3]));
        mxp(~sum(mxp,2),:) = 1/3;  %correct rows that sum zero
        mxp = diag(1./sum(mxp,2))*mxp;   %re-adjusting to sum one
        mxp(:,4) = 0;

        %%% INSERT TRANSITION PROB
        ixp = exp(-XState(:,[4 5]));
        ixp(~sum(ixp,2),:) = 1/2;  %correct rows that sum zero
        ixp = diag(1./sum(ixp,2))*ixp;   %re-adjusting to sum one

        %%% DELETE TRANSITION PROB
        dxp = exp(-XState(:,[6 7]));
        dxp(~sum(dxp,2),:) = 1/2;  %correct rows that sum zero
        dxp = diag(1./sum(dxp,2))*dxp;   %re-adjusting to sum one

        %%% LEFT FLANKING INSERT TRANSITION PROB
        lixp = nxp;

        %%% RIGTH FLANKING INSERT TRANSITION PROB
        rixp = nxp;
        
        %%% LOOP INSERT TRANSITION PROB
        lxp = nxp;
        
        %%% END TRANSITION PROB
        %%% [E->C, E->J]
        endxp = [.5 .5];

        %%% Store data in structure
        data.MatchEmission = mep;
        data.InsertEmission = iep;
        data.NullEmission = nep;
        data.BeginX = bxp';
        data.MatchX = mxp(1:end-1,:); %last state is always [0 0 0 1]
        data.InsertX = ixp(1:end-1,:); %last state does not exist
        data.DeleteX = dxp(1:end-1,:); %last state always goes to E
        data.FlankingInsertX = [lixp;rixp]';
        data.LoopX = [endxp;lxp]';
        data.NullX = nxp';
end

function str = removeallblanks(str)
% REMOVEBLANKS removes all blanks
str(isspace(str))=[];
end

function out = CorrectStars(in)
% CORRECTSTARS correct the string marked with '*' which represents Zero
% probabilities (i.e. -Inf log-Probabilities), returns the numeric array
out = sscanf(strrep(in,'*','-Inf'),'%f');
end

function out = CorrectStars3b(in)
% CORRECTSTARS correct the string marked with '*' which represents Zero
% probabilities (i.e. -Inf log-Probabilities), returns the numeric array
out = sscanf(strrep(in,'*','Inf'),'%f');
end
