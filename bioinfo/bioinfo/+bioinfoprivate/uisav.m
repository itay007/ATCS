classdef uisav < hgsetget
    % version v0.94
    % Author: elmar.tarajan@mathworks.de
    
    properties ( GetAccess = protected )
        % Only properties that can be set by the client,
        % the client can not depend on the state of these properties,
        % therefore there is not read access.
        Data
        Scheme
        FontSize = getArchDependentDefaultFontSize
        EnableBackground = 1
        
    end
    properties (Dependent, GetAccess = protected )
        ContextMenu
    end
    
    properties (Constant, Access = protected )
        MouseIcons = load(fullfile(toolboxdir('bioinfo'),'bioinfo','+bioinfoprivate','seqAlignViewerIcons.mat'),'OpenHand','CloseHand')
        FontName = getArchDependentFontName
        FontProps = getArchDependentFontProperties(getArchDependentFontName)
    end
    
  
    properties ( Access = protected )

        Parent
        hFig
        
        hMainPanel
        
        hAxesConsensus
        hTextConsensus
        hBarPlotConsensus
        
        hAxesAlignment
        hTextAlignment
        hBackgroundAlignment
        hSelectBox
        hSelectRows
        hSelectColumns
        
        hAxesMagnifier
        hBackgroundMagnifier
        hBoxMagnifier
        hBoxTopLeftText
        hBoxTopRightText
        hBoxLeftTopText
        hBoxLeftBottomText
        
        hSliderHorizontal
        hSliderVertical
        hFooterRight
        hFooterLeft
        
        seqWin
        seqInd
        cellSize = [10 18]
        winSize
        xWin = 0
        yWin = 0
        x
        y

        xSel
        ySel

        xInd
        yInd
        SeqLabelExtent = 100
        
        seq_txt
        seq_names
        seq_rgb
        seq_consensus
        seq_scores
        
        undo
        redo
        
        listener
    end
    
    
    methods
        function obj = uisav(varargin)
            % parse input parameter
            par = inputParser;
            par.addParamValue('Parent',[]);
            par.addParamValue('Data',[]);
            par.addParamValue('Scheme',struct('Name','Gray','Symbols','ABCDEFGHIKLMNPQRSTUVWYZX','Colors',gray(24)*0.5+.5));
            par.parse(varargin{:})
            par = par.Results;
            
            % build-up the user interface
            if isempty(par.Parent)
                par.Parent = figure;
            end% if
            obj.hMainPanel = uipanel( ...
                'Tag','uisav.MainPanel',...
                'Parent',par.Parent, ...
                'units','normalized', ...
                'bordertype','none', ...
                'position',[0 0 1 1], ...
                'bordertype','none', ...
                'Highlightcolor',[.5 .5 .5], ...
                'visible','on');
            obj.hFig = ancestor(obj.hMainPanel,'figure');
            obj.Parent = par.Parent;
            
            % consensus axes
            obj.hAxesConsensus = axes( ...
                'Tag','uisav.AxesConsensus',...
                'Parent',obj.hMainPanel, ...
                'units','pixel', ...
                'nextplot','add', ...
                'FontSize',7, ...
                'fontunits','pixel', ...
                'FontWeight','bold', ...
                'layer','top', ...
                'box','on', ...
                'drawmode','fast', ...
                'color',[ 1.0000    0.8627    0.5882 ], ...
                'XColor',[.2 .2 .2], ...
                'YColor',[.2 .2 .2], ...
                'XTickLabel','', ...
                'YTickLabel',['' 'Consensus'], ...
                'TickLength',[0 0], ...
                'XAxisLocation','bottom', ...
                'ButtonDownFcn',{@obj.MouseAction,'down_in_consensus_axes'});
            
            obj.hTextConsensus = text(0,0,{'  ' '  '}, ...
                'Tag','uisav.TextConsensus',...
                'Parent',obj.hAxesConsensus, ...
                'verticalalignment','bottom', ...
                'horizontalalignment','left', ...
                'units','pixel', ...
                'fontunits','pixel', ...
                'fontname',obj.FontName, ...
                'units','pixel', ...
                'FontSize',obj.FontSize, ...
                'FontWeight','normal', ...
                'HitTest','off', ...
                'clipping','on');
            
            % consensus barplot
            obj.hBarPlotConsensus = patch(NaN,NaN,[0 0 1], ...
                'Tag','uisav.BarPlotConsensus',...
                'Parent',obj.hAxesConsensus, ...
                'edgeColor','none', ...
                'FaceColor',[0.7843 0.3922 0.1961]);
            
            % main axes
            obj.hAxesAlignment = axes( ...
                'Tag','uisav.AxesAlignment',...
                'Parent',obj.hMainPanel, ...
                'units','pixel', ...
                'nextplot','add', ...
                'FontSize',obj.FontSize, ...
                'fontunits','pixel', ...
                'FontWeight','normal', ...
                'TickDir','out', ...
                'YDir','reverse', ...
                'layer','top', ...
                'box','on', ...
                'drawmode','fast', ...
                'color','none', ...
                'XColor',[.2 .2 .2], ...
                'YColor',[.2 .2 .2], ...
                'TickLength',[0 0], ...
                'XAxisLocation','bottom', ...
                'ButtonDownFcn',{@obj.MouseAction,'down_in_main_axes_labels'});
            
            % backgroundcolor-layer
            obj.hBackgroundAlignment = image(nan(1), ...
                'Tag','uisav.BackgroundAlignment',...
                'Parent',obj.hAxesAlignment, ...
                'Uicontextmenu',[], ...
                'buttondownfcn',{@obj.MouseAction,'down_inside_main_axes'});
            
            % text-layer
            obj.hTextAlignment = text(0,0,{'  ' '  '}, ...
                'Tag','uisav.TextAlignment',...
                'Parent',obj.hAxesAlignment, ...
                'verticalalignment','bottom', ...
                'horizontalalignment','left', ...
                'units','pixel', ...
                'fontunits','pixel', ...
                'fontname',obj.FontName, ...
                'units','pixel', ...
                'FontSize',obj.FontSize, ...
                'FontWeight','normal', ...
                'HitTest','off', ...
                'clipping','on');
            
            % selection patches
            obj.hSelectRows = patch(NaN,NaN,[0 0 1], ...
                'Tag','uisav.SelectRows',...
                'Parent',obj.hAxesAlignment, ...
                'facecolor','none', ...
                'EdgeColor',[.4 .4 .4], ...
                'HitTest','on', ...
                'lineWidth',1,...
                'buttondownfcn',{@obj.MouseAction,'down_in_selection'});
            obj.hSelectColumns = patch(NaN,NaN,[0 0 1], ...
                'Tag','uisav.SelectColumns',...
                'Parent',obj.hAxesAlignment, ...
                'facecolor','none', ...
                'EdgeColor',[.4 .4 .4], ...
                'HitTest','off', ...
                'lineWidth',1);
            obj.hSelectBox = patch(NaN,NaN,NaN, ...
                'Tag','uisav.SelectBox',...
                'Parent',obj.hAxesAlignment, ...
                'EdgeColor',[1 0 0], ...
                'HitTest','off', ...
                'Uicontextmenu',[], ...
                'LineWidth',1);
            % Axes for whole data and magnifier-patch
            obj.hAxesMagnifier = axes( ...
                'Tag','uisav.AxesMagnifier',...
                'Parent',obj.hMainPanel, ...
                'units','pixel', ...
                'YDir','reverse', ...
                'layer','top', ...
                'XTickLabel',[], ...
                'YTickLabel',[], ...
                'TickLength',[0 0], ...
                'box','on', ...
                'visible','off', ...
                'nextplot','add');
            obj.hBackgroundMagnifier = image(NaN,...
                'Tag','uisav.BackgroundMagnifier',...
                'Parent',obj.hAxesMagnifier); % ,'buttondownfcn',@SetWindow);
            obj.hBoxMagnifier = patch(NaN,NaN,[1 0 0], ...
                'Tag','uisav.BoxMagnifier',...
                'Parent',obj.hAxesMagnifier, ...
                'FaceColor','none', ...
                'Edgecolor',[1 0 0], ...
                'LineWidth',1, ...
                'buttondownfcn',{@obj.MouseAction,'down_in_magnifier'});
            obj.hBoxTopLeftText = text(NaN,NaN,'', ...
                'Parent',obj.hAxesMagnifier, ...
                'Color',[0 0 0], ...
                'fontunits','pixel', ...
                'fontname',obj.FontName, ...
                'fontsize',obj.FontSize-1, ...
                'FontWeight','normal', ...
                'HorizontalAlignment','left', ...
                'VerticalAlignment','bottom', ...
                'HitTest','off');
            obj.hBoxTopRightText = text(NaN,NaN,'', ...
                'Parent',obj.hAxesMagnifier, ...
                'Color',[0 0 0], ...
                'fontunits','pixel', ...
                'fontname',obj.FontName, ...
                'fontsize',obj.FontSize-1, ...
                'FontWeight','normal', ...
                'HorizontalAlignment','right', ...
                'VerticalAlignment','bottom', ...
                'HitTest','off');
            obj.hBoxLeftTopText = text(NaN,NaN,'', ...
                'Parent',obj.hAxesMagnifier, ...
                'Color',[0 0 0], ...
                'fontunits','pixel', ...
                'fontname',obj.FontName, ...
                'fontsize',obj.FontSize-1, ...
                'FontWeight','normal', ...
                'HorizontalAlignment','right', ...
                'VerticalAlignment','top', ...
                'HitTest','off');
            obj.hBoxLeftBottomText = text(NaN,NaN,'', ...
                'Parent',obj.hAxesMagnifier, ...
                'Color',[0 0 0], ...
                'fontunits','pixel', ...
                'fontname',obj.FontName, ...
                'fontsize',obj.FontSize-1, ...
                'FontWeight','normal', ...
                'HorizontalAlignment','right', ...
                'VerticalAlignment','bottom', ...
                'HitTest','off');
            obj.hFooterLeft = uicontrol( ...
                'Tag','uisav.FooterLeft',...
                'Parent',obj.hMainPanel, ...
                'style','checkbox', ...
                'cdata',nan(1,1,3), ...
                'units','pixel', ...
                'enable','inactive', ...
                'FontSize',obj.FontSize-1, ...
                'FontName',obj.FontName, ...
                'HorizontalAlignment','left', ...
                'foregroundcolor',[.4 .4 .4]);
            obj.hFooterRight = uicontrol( ...
                'Tag','uisav.FooterRight',...
                'Parent',obj.hMainPanel, ...
                'style','edit', ...
                'units','pixel', ...
                'enable','inactive', ...
                'FontSize',obj.FontSize-1, ...
                'FontName',obj.FontName, ...
                'HorizontalAlignment','right', ...
                'foregroundcolor',[.4 .4 .4]);            
            obj.hSliderHorizontal = uicontrol( ...
                'Tag','uisav.SliderHorizontal',...
                'Parent',obj.hMainPanel, ...
                'style','slider', ...
                'units','pixel', ...
                'min',1, ...
                'max',2, ...
                'value',1);
            obj.hSliderVertical = uicontrol( ...
                'Tag','uisav.SliderVertical',...
                'Parent',obj.hMainPanel, ...
                'style','slider', ...
                'units','pixel', ...
                'min',1, ...
                'max',2, ...
                'value',1);            
            try % for newer matlab versions
                addlistener(obj.hSliderHorizontal,'ContinuousValueChange',@obj.UpdateContent);
                addlistener(obj.hSliderVertical,'ContinuousValueChange',@obj.UpdateContent);
            catch ME %#ok<NASGU>  for older matlab verisons
                addlistener(obj.hSliderHorizontal,'Action',@obj.UpdateContent);
                addlistener(obj.hSliderVertical,'Action',@obj.UpdateContent);
            end
            
            
            % set Scheme
            obj.Scheme = par.Scheme;
            
            % set data
            obj.Data = par.Data;
            
            % Keyboard support (TBD)
            % ----------------------
            %   set(obj.hFig,'KeyPressFcn',@obj.KeyPress)
            %   set(obj.hFig,'KeyReleaseFcn',@obj.KeyRelease)
            
            % initialize listener
            % -------------------
            obj.listener{1} = addlistener(obj.hMainPanel,'PixelBounds','PostSet',@obj.UpdatePositions);

            obj.FontSize = getArchDependentDefaultFontSize;            
            obj.UpdateContent
            obj.FontSize = getArchDependentDefaultFontSize;  
        end
        
        function set.ContextMenu(obj,value)
            set(obj.hSelectRows,'Uicontextmenu',value)
        end
        
        function set.EnableBackground(obj,value)
            switch value
                case {'on' 1}
                    value = 1;
                case {'off' 0}
                    value = 0;
            end% switch
            obj.EnableBackground = value;
            if ~obj.IsDataEmpty
                obj.Select(obj.xSel,obj.ySel);
            end
        end
        
        function set.Scheme(obj,value)
            if ~isa(value,'struct')
                fprintf(2,'Input value is not a valid color scheme structure\n')
            else
                if ~all(ismember(fieldnames(value)',{'Name' 'Type' 'Symbols' 'Colors'}))
                    fprintf(2,'Input value is not a valid color scheme structure\n')
                else
                    obj.Scheme = value;
                    set(obj.hFooterLeft,'String',obj.Scheme.Name)
                    ext = get(obj.hFooterLeft,'Extent');
                    pos = get(obj.hFooterLeft,'Position');
                    pos(3) = ext(3)+8;
                    set(obj.hFooterLeft,'position',pos)
                    
                    clr_label = ['.' obj.Scheme.Symbols];
                    clr_value = [1 1 1 ; obj.Scheme.Colors];
                    seq_clr = zeros(size(obj.seq_txt));
                    for n = 1:numel(clr_label)
                        seq_clr(obj.seq_txt==clr_label(n)) = n;
                    end% for
                    %  seq_clr(seq_clr==0) = n+1;
                    obj.seq_rgb = ind2rgb(seq_clr,clr_value);
                    
                    % update magnifier axes
                    set(obj.hBackgroundMagnifier,'CData',obj.seq_rgb)
                    
                    if ~obj.IsDataEmpty
                        obj.Select(obj.xSel,obj.ySel);
                    end
                end
            end
        end
        
        function set.FontSize(obj,value)
            if isempty(value)
                return
            else
                obj.FontSize = obj.FontProps.FontSizeReel(find([(value<=obj.FontProps.FontSizeReel(1:end-1)) true],1,'first'));
            end
                     
            % calculate character size in pixel
            %set(obj.hTextAlignment,'FontSize',obj.FontSize,'String',repmat(' ',2,2))
            set(obj.hTextAlignment,'FontSize',obj.FontSize)
            tmp = get(obj.hTextAlignment,'extent');
            obj.cellSize = (tmp(3:4)-obj.FontProps.ExtentEdge)./fliplr(size(get(obj.hTextAlignment,'String'))).*2;

            if isempty(get(obj.hTextAlignment,'String'))
                set(obj.hTextAlignment,'FontSize',obj.FontSize,'String',repmat(' ',2,2));
                tmp = get(obj.hTextAlignment,'extent');
                obj.cellSize = tmp(3:4)-obj.FontProps.ExtentEdge;
            end
            
            % set new fontsize
            set([obj.hTextAlignment obj.hTextConsensus obj.hAxesAlignment obj.hAxesConsensus],'FontSize',obj.FontSize)
            set([obj.hTextAlignment obj.hTextConsensus],'position',[0 2]+(obj.cellSize./[4 -4]));
            
            % determine extent for sequence names
            set(obj.hTextAlignment,'FontSize',obj.FontSize,'String',char(obj.seq_names))
            obj.SeqLabelExtent = get(obj.hTextAlignment,'extent');
            obj.SeqLabelExtent = obj.SeqLabelExtent(3)+25;
            
            % update
            if ~obj.IsDataEmpty
                obj.UpdatePositions
            end
        end
        
        function set.Data(obj,value)
            if ~isempty(value)
                obj.Data = value;
                % create data arrays
                obj.seq_txt = upper(char({obj.Data.Sequence}));
                obj.seq_txt(obj.seq_txt==' ')='.';
                obj.seq_txt(obj.seq_txt=='-')='.'; % Most common symbol used to represent gaps
                obj.seq_txt(obj.seq_txt=='_')='.';
                obj.seq_names = {obj.Data.Header};
                obj.x = size(obj.seq_txt,2);
                obj.y = size(obj.seq_txt,1);
                
                clr_label = ['.' obj.Scheme.Symbols];
                clr_value = [1 1 1 ; obj.Scheme.Colors];
                seq_clr = zeros(size(obj.seq_txt));
                for n = 1:numel(clr_label)
                    seq_clr(obj.seq_txt==clr_label(n)) = n;
                end% for
                %  seq_clr(seq_clr==0) = n+1;
                obj.seq_rgb = ind2rgb(seq_clr,clr_value);
                
                % update magnifier axes
                set(obj.hBackgroundMagnifier,'CData',obj.seq_rgb)
                set(obj.hAxesMagnifier,'xlim',[0 obj.x]+0.5,'ylim',[0 obj.y]+0.5)
                
                % determine extent for sequence names
                set(obj.hTextAlignment,'FontSize',obj.FontSize,'String',char(obj.seq_names))
                obj.SeqLabelExtent = get(obj.hTextAlignment,'extent');
                obj.SeqLabelExtent = obj.SeqLabelExtent(3)+30;
                
                % update
                obj.UpdatePositions
                
                % remove selection
                obj.Select([],[]);
                
                % init mouse
                set(obj.hFig, ...
                    'WindowButtonMotionFcn',{@obj.MouseAction,'mouse_hover'}, ...
                    'WindowButtonUpFcn',{@obj.MouseAction,'up'})
                
            end
        end
        
        function Select(obj,xPos,yPos)
            % Select/Highlight area
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Select(5:10,1:10) % Select block
            if nargin==2
                if xPos % Select(1) ~ Select all, especial signature used by the client
                    obj.Select(1:obj.x,1:obj.y)
                else    % Select(0) ~ Deselect all, especial signature used by the client
                    obj.Select([],[]);
                end
                return
            end% if
            % there can never exist a selection in one direction while the
            % other direction is empty
            if isempty(xPos) && ~isempty(yPos)
                xPos = 1:obj.x;
            elseif ~isempty(xPos) && isempty(yPos)
                yPos = 1:obj.y;
            end
            
            % check limits
            obj.xSel = max(1,min(xPos,obj.x));
            obj.ySel = max(1,min(yPos,obj.y));
            %
            % highlight area in the magnifier axes
            rgb = obj.seq_rgb;
            rgb(obj.ySel,obj.xSel,:) = rgb(obj.ySel,obj.xSel,:)*.75;
            set(obj.hBackgroundMagnifier,'CData',rgb);
            
            % highlight main view
            switch obj.EnableBackground
                case 0
                    % white rgb array
                    rgb = ones(numel(obj.yInd),numel(obj.xInd),3);
                    rgb(ismember(obj.yInd,obj.ySel),ismember(obj.xInd,obj.xSel),:) = 0.75;
                    set(obj.hBackgroundAlignment,'CData',rgb)
                case 1
                    set(obj.hBackgroundAlignment,'CData',rgb(obj.yInd,obj.xInd,:))
            end% switch
            obj.UpdateSelection
        end
        
        function View(obj,y,x)
            % Set magnifier rectangle (left upper corner)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % calculate new values for horizontal and vertical sliders
            value = min(max(x,get(obj.hSliderHorizontal,'min')),get(obj.hSliderHorizontal,'max'));
            set(obj.hSliderHorizontal,'value',value);
            % set vertical slider
            value = get(obj.hSliderVertical,'max')-min(max(y,get(obj.hSliderVertical,'min')),get(obj.hSliderVertical,'max'))+1;
            set(obj.hSliderVertical,'value',value);
            %
            obj.UpdateContent(obj)
        end
        
        function Order(obj,yInd)
            % Reorder the sequences
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if all(ismember(1:numel(obj.Data),yInd))
                obj.Data = obj.Data(yInd);
            end
        end
        
        function tf = IsDataEmpty(obj)
            % Returns true/false to indicate if there is not data contained
            % in the uisav (used by the client to update menus on the fly)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tf = isempty(obj.Data);
        end
        
        function tf = IsAreaSelected(obj)
            % Returns true/false to indicate if there is an active selected
            % area (used by the client to update menus on the fly)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tf = numel(obj.ySel)>0;
        end
        
        function fs = getFontSize(obj)
            % Returns the current font size
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fs = obj.FontSize;
        end
        function fsr = getFontSizeReel(obj)
            % Returns the current font size
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fsr = obj.FontProps.FontSizeReel;
        end        
        
        function output = ExportSelectedArea(obj,yPos,xPos)
            % get selected area as struct
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin==1
                yPos = obj.ySel;
                xPos = obj.xSel;
            end% if
            if isempty(yPos) % if nothing selected the whole data is exported
                yPos = 1:obj.y;
                xPos = 1:obj.x;
            end
            output = repmat(struct('Header',[],'Sequence',[]),numel(yPos),1);
            for nRow = 1:numel(yPos)
                output(nRow,:).Header   = obj.seq_names{yPos(nRow)};
                output(nRow,:).Sequence = strrep(obj.seq_txt(yPos(nRow),xPos),'.','-');
            end
        end
        
        function output = ExportConsensus(obj,xPos)
            % get consensus as a structure
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin==1
                xPos = obj.xSel;
            end
            if isempty(xPos) % if nothing selected the whole consensus is exported
                ranges = [1,obj.x];
            else % discontinuos selections are exported as multiple sequences
                xd = diff([0;accumarray(xPos(:),1);0]);
                ranges = [find(xd>0) find(xd<0)-1];
            end
            output = repmat(struct('Header',[],'Sequence',[]),size(ranges,1),1);
            for nRow = 1:size(ranges,1)
                if diff(ranges(nRow,:))
                    output(nRow,:).Header = sprintf('Consensus %d-%d',ranges(nRow,:));
                else
                    output(nRow,:).Header = sprintf('Consensus %d',ranges(nRow,1));
                end
                output(nRow,:).Sequence = seqconsensus(obj.seq_txt(:,ranges(nRow,1):ranges(nRow,2)),'GAPS','ALL');
            end
        end
        
        
        function DeleteSequences(obj,rows)
            % Delete whole rows of the selection
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin==1
                rows = obj.ySel;
                if ~isempty(rows)
                    answer = questdlg(sprintf('Delete the selected %d Sequences?',numel(rows)), ...
                        'Delete','Yes', 'No', 'No');
                    switch answer
                        case 'Yes'
                            rows = obj.ySel;
                        case 'No'
                            rows = [];
                    end % switch
                end % if
            end
            if ~isempty(rows)
                obj.Data(rows) = [];
                obj.ySel = [];
                obj.xSel = [];
            end
        end
        
        function DeleteColumns(obj,cols)
            % Delete selected columns that are all gaps
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin==1
                cols = obj.xSel;
            end
            if isempty(cols) % if nothing selected the whole data is checked for empty cols
                cols = 1:obj.x;
            end
            
            tmp = ExportSelectedArea(obj,[],[]); % gets all the data
            cols = cols(all(subsref(char(tmp.Sequence),substruct('()',{':',cols}))=='-',1));
            if ~isempty(cols)
                for nSeq = 1:obj.y
                    tmp(nSeq).Sequence(cols) = [];
                end
                xSeltmp = obj.xSel;
                ySeltmp = obj.ySel;
                % Updates the data but keeps the selection, shrink the
                % selection as necessary...
                obj.Data = tmp;
                if ~isempty(xSeltmp)
                    acc = cumsum(accumarray(cols(:),1,[max(xSeltmp),1]))';
                    obj.xSel = setdiff(xSeltmp,cols)-acc(setdiff(xSeltmp,cols));
                    if ~isempty(obj.xSel)
                        obj.ySel = ySeltmp;
                    end
                end
            end
        end
        
        
        function varargout = Undo(obj)
            if nargout
                varargout{1} = numel(obj.undo);
                varargout{2} = obj.undo;
            else
                if ~isempty(obj.undo)
                    obj.Select(obj.undo{end}{2},obj.undo{end}{1})
                    [ xpos,step ] = obj.MoveBlock(obj.undo{end}{:});
                    switch obj.undo{end}{4}
                        case 1
                            obj.redo{end+1} = {obj.undo{end}{1} xpos -step obj.undo{end}{4}};
                            obj.Select(xpos,obj.undo{end}{1})
                        case 2
                            obj.redo{end+1} = {xpos obj.undo{end}{2} -step obj.undo{end}{4}};
                            obj.Select(obj.undo{end}{2},xpos)
                    end
                    obj.undo(end) = [];
                    %
                    obj.UpdateContent
                    %
                    drawnow
                end
            end
        end
        
        function varargout = Redo(obj)
            if nargout
                varargout{1} = numel(obj.redo);
                varargout{2} = obj.redo;
            else
                if ~isempty(obj.redo)
                    [ xpos,step ] = obj.MoveBlock(obj.redo{end}{:});
                    switch obj.redo{end}{4}
                        case 1
                            obj.undo{end+1} = {obj.redo{end}{1} xpos -step obj.redo{end}{4}};
                            obj.Select(xpos,obj.redo{end}{1})
                        case 2
                            obj.undo{end+1} = {xpos obj.redo{end}{2} -step obj.redo{end}{4}};
                            obj.Select(obj.redo{end}{2},xpos)
                    end
                    obj.redo(end) = [];
                    obj.UpdateContent
                    drawnow
                end
            end
        end
        
        function MoveSelectedBlock(obj,step,direction)
            % Moves the current selected block (this method is used by the
            % client; as the client should rather not have access to obj.xSel
            % and obj.ySel)
            % step may be -inf or -inf to move all the way to the limits or
            % any other integer.
            % direction must be 1 or 2.
            
            switch direction
                case 1 % horz
                    step = min(max(step,1-min(obj.xSel)),obj.x-max(obj.xSel));
                    obj.MoveBlock(obj.ySel,obj.xSel,step,1);
                case 2 % vert
                    step = min(max(step,1-min(obj.ySel)),obj.y-max(obj.ySel));
                    obj.MoveBlock(obj.ySel,obj.xSel,step,2);
            end
        end
        
        function [fpos,step] = MoveBlock(obj,ypos,xpos,step,direction)
            switch direction % horz
                case 1
                    if isequal(xpos,obj.xSel)
                        updateXSel = true;
                    else
                        updateXSel = false;
                    end
                    blk = obj.seq_txt(ypos,xpos);
                    clr = obj.seq_rgb(ypos,xpos,:);
                    if step > 0
                        pos = xpos(end)+(1:abs(step));
                        tmp = all(obj.seq_txt(ypos,pos(pos<=size(obj.seq_txt,2)))=='.',1);
                        if all(tmp)
                            addlines = sum(pos>size(obj.seq_txt,2));
                            obj.seq_txt = [obj.seq_txt repmat('.',size(obj.seq_txt,1),addlines)];
                            obj.seq_rgb = [obj.seq_rgb ones(size(obj.seq_txt,1),addlines,3)];
                        else
                            step = find(tmp==0,1,'first')-1;
                        end
                        obj.seq_txt(ypos,(xpos(1)):(xpos(end)+step)) = [repmat('.',numel(ypos),abs(step)) blk];
                        obj.seq_rgb(ypos,(xpos(1)):(xpos(end)+step),:) = [ones(numel(ypos),abs(step),3) clr];
                        xpos = xpos + step;
                        
                        % delete gaps if necessary
                        ind = xpos(1)-1;
                        if ind
                            ind = find([all(obj.seq_txt(:,1:ind)=='.',1) 0]==0,1,'first')-1;
                            obj.seq_txt(:,1:ind)   = [];
                            obj.seq_rgb(:,1:ind,:) = [];
                            xpos = xpos - ind;
                            obj.xSel = obj.xSel - ind;
                        end
                        
                    elseif step < 0
                        pos = xpos(1)-(1:abs(step));
                        tmp = all(obj.seq_txt(ypos,pos(pos>0))=='.',1);
                        if all(tmp)
                            addlines = sum(pos<=0);
                            obj.seq_txt = [repmat('.',size(obj.seq_txt,1),addlines) obj.seq_txt];
                            obj.seq_rgb = [ones(size(obj.seq_txt,1),addlines,3) obj.seq_rgb];
                        else
                            addlines = 0;
                            step = -find(tmp==0,1,'first')+1;
                        end
                        obj.seq_txt(ypos,(xpos(1)+addlines+step):(xpos(end)+addlines)) = [blk repmat('.',numel(ypos),abs(step))];
                        obj.seq_rgb(ypos,(xpos(1)+addlines+step):(xpos(end)+addlines),:) = [clr ones(numel(ypos),abs(step),3)];
                        xpos = xpos + addlines + step;
                        obj.xSel = obj.xSel + addlines ;
                        
                        % delete gaps if necessary
                        ind = size(obj.seq_txt,2)-xpos(end);
                        if ind
                            ind = find([all(obj.seq_txt(:,end:-1:end-ind+1)=='.',1) 0]==0,1,'first')-1;
                            obj.seq_txt(:,end:-1:end-ind+1)   = [];
                            obj.seq_rgb(:,end:-1:end-ind+1,:) = [];
                        end
                    end
                    
                    if step~=0
                        %obj.xSel = xpos;
                        % update uisav-object
                        [obj.y,obj.x] = size(obj.seq_txt);
                        set(obj.hBackgroundMagnifier,'CData',obj.seq_rgb);
                        set(obj.hAxesMagnifier,'xlim',[0 obj.x]+0.5,'ylim',[0 obj.y]+0.5)
                    end
                    
                    fpos = xpos;
                    if updateXSel
                        obj.xSel = xpos;
                    end
                    
                    obj.UpdatePositions
                    
                    obj.Select(obj.xSel,obj.ySel)
                case 2 % vertical
                    step = min(max(step,1-min(ypos)),obj.y-max(ypos));
                    % find permuting vector
                    perm = zeros(obj.y,1);
                    perm(ypos+step) = ypos;
                    perm(perm==0) = setdiff(1:obj.y,ypos);
                    % update data
                    obj.seq_txt = obj.seq_txt(perm,:);
                    obj.seq_rgb = obj.seq_rgb(perm,:,:);
                    obj.seq_names = obj.seq_names(perm);
                    set(obj.hBackgroundMagnifier,'CData',obj.seq_rgb);
                    set(obj.hBackgroundAlignment,'CData',obj.seq_rgb(obj.yInd,obj.xInd,:))
                    obj.seqWin(obj.seqInd) = obj.seq_txt(obj.yInd,obj.xInd);
                    set(obj.hTextAlignment,'String',obj.seqWin)
                    obj.UpdatePositions % rendering update
                    % new positions of rows being moved
                    fpos = ypos + step; % put resulting pos in the output argument
                    obj.Select(obj.xSel,fpos) % update highligthed selection
            end %switch direction
        end % function MoveBlock
    end % methods
    
    methods (Access = protected)
        
        function MouseAction(obj,~,~,action)
            % Mouse control for all interactive activities
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            persistent posG0 offset xblock
            % Selection Type (Alt,Extend,Normal,Open)
            selType = get(obj.hFig,'selectionType');
            % mouse position
            pos = round(get(obj.hAxesAlignment,'CurrentPoint'));
            % mouse global position
            posG = [obj.xInd(1) obj.yInd(1)]+pos([1 3])-1;
            posG(1) = min(obj.x,max(1,posG(1)));
            posG(2) = min(obj.y,max(1,posG(2)));
            % determine if mouse is within the axes (true, false)
            ptrInAxes = all([pos([1 3])>0 [obj.xWin obj.yWin]-pos([1 3])>=0]);
            % update footer
            obj.UpdateFooter(ptrInAxes,posG)
            
            % debug me ....
            % fprintf('Action: %s \t selType: %s \n',action,selType)
            
            %
            switch action
                case 'down_inside_main_axes'
                    switch selType
                        case 'normal'
                            posG0 = posG;
                            set(obj.hFig,'WindowButtonMotionFcn',{@obj.MouseAction,'selecting_region_unrestricted'})
                            obj.Select(posG0(1),posG0(2))
                        case 'alt'
                            obj.Select([],[]) %resets selection
                    end
                    
                case 'down_in_selection'
                    switch selType
                        case 'normal' % start to move the block left-right
                            set(obj.hFig,'PointerShapeCData',obj.MouseIcons.CloseHand)
                            posG0 = posG;
                            if posG(1)<min(obj.xSel)
                                xblock = 1:(min(obj.xSel)-1);
                            elseif posG(1)>max(obj.xSel)
                                xblock = (max(obj.xSel)+1):obj.x;
                            else
                                xblock = min(obj.xSel):max(obj.xSel);
                                if ~all(diff(obj.xSel)==1)
                                    return
                                end
                            end
                            if ~isempty(xblock)
                                set(obj.hFig, ...
                                    'WindowButtonMotionFcn',{@obj.MouseAction,'moving_block_left_right'}, ...
                                    'WindowButtonUpFcn',{@obj.MouseAction,'up'})
                            end
                            obj.MouseAction([],[],'moving_block_left_right')
                            
                        case 'extend' % start to move the block up-down
                            set(obj.hFig,'PointerShapeCData',obj.MouseIcons.CloseHand)
                            posG0 = posG;
                            set(obj.hFig,'WindowButtonMotionFcn',{@obj.MouseAction,'moving_block_up_down'})
                            
                        case 'alt' % reset selection
                            % reset selection by right mouse click outside the selected area
                            if ~(any(obj.ySel==posG(2)) && any(obj.xSel==posG(1)))
                                obj.Select([],[]) %resets selection
                            end
                    end
                    
                case 'down_in_consensus_axes'
                    if any(strcmp(selType,{'alt' 'open'}))
                        xPos = setxor(obj.xSel,posG(1));
                        obj.Select(xPos,obj.ySel);
                    else
                        posG0 = posG;
                        set(obj.hFig,'WindowButtonMotionFcn',{@obj.MouseAction,'selecting_region_y_restricted'})
                        obj.MouseAction([],[],'selecting_region_y_restricted')
                    end
                    
                case 'down_in_main_axes_labels'
                    if any(strcmp(selType,{'alt' 'open'}))
                        yPos = setxor(obj.ySel,posG(2));
                        obj.Select(obj.xSel,yPos);
                    else                         posG0 = posG;
                        set(obj.hFig,'WindowButtonMotionFcn',{@obj.MouseAction,'selecting_region_x_restricted'})
                        obj.MouseAction([],[],'selecting_region_x_restricted')
                    end
                    
                case 'selecting_region_x_restricted'
                    switch selType
                        case 'normal'
                            if isempty(obj.xSel)
                                obj.xSel = 1:obj.x;
                            end
                            obj.ySel = [];
                            yPos = union(obj.ySel,min(posG(2),posG0(2)):max(posG(2),posG0(2)));
                        case 'extend'
                            yPos = union(obj.ySel,min(posG(2),posG0(2)):max(posG(2),posG0(2)));
                            %                      case 'alt'
                            %                         yPos = setxor(obj.ySel,posG(2));
                    end
                    obj.Select(obj.xSel,yPos)
                    
                case 'selecting_region_y_restricted'
                    switch selType
                        case 'normal'
                            if isempty(obj.ySel)
                                obj.ySel = 1:obj.y;
                            end
                            obj.xSel = [];
                            xPos = union(obj.xSel,min(posG(1),posG0(1)):max(posG(1),posG0(1)));
                        case 'extend'
                            xPos = union(obj.xSel,min(posG(1),posG0(1)):max(posG(1),posG0(1)));
                            %                      case 'alt'
                            %                         xPos = setxor(obj.xSel,posG(1));
                    end
                    obj.Select(xPos,obj.ySel)
                    
                case 'moving_block_left_right'
                    if ~isequal(posG0(1),posG(1))
                        xblock = obj.MoveBlock(obj.ySel,xblock,posG(1)-posG0(1),1);
                        posG0 = posG;
                    end
                    
                case 'moving_block_up_down'
                    if ~isequal(posG0(2),posG(2))
                        obj.MoveBlock(obj.ySel,obj.xSel,posG(2)-posG0(2),2);
                        posG0 = posG;
                    end
                    
                case 'selecting_region_unrestricted'
                    xPos = min(posG(1),posG0(1)):max(posG(1),posG0(1));
                    yPos = min(posG(2),posG0(2)):max(posG(2),posG0(2));
                    obj.Select(xPos,yPos)
                    
                case 'mouse_hover'
                    if ptrInAxes
                        if any(obj.ySel==posG(2))% && any(obj.xSel==posG(1))
                            set(obj.hFig,'Pointer','Custom','PointerShapeCData',obj.MouseIcons.OpenHand)
                        else
                            set(obj.hFig,'Pointer','Arrow')
                        end
                    else
                        pos = round(get(obj.hAxesMagnifier,'CurrentPoint'));
                        xdata = get(obj.hBoxMagnifier,'xdata');
                        ydata = get(obj.hBoxMagnifier,'ydata');
                        if (pos(1)>min(xdata) && pos(1)<max(xdata)) && (pos(3)>min(ydata) && pos(3)<max(ydata))
                            set(obj.hFig,'Pointer','Custom','PointerShapeCData',obj.MouseIcons.OpenHand)
                        else
                            set(obj.hFig,'Pointer','Arrow')
                        end
                    end
                    
                case 'down_in_magnifier'
                    set(obj.hFig,'PointerShapeCData',obj.MouseIcons.CloseHand)
                    posM = round(get(obj.hAxesMagnifier,'CurrentPoint'));
                    ydata = get(obj.hBoxMagnifier,'ydata');
                    xdata = get(obj.hBoxMagnifier,'xdata');
                    offset = round(posM([1 3])-[xdata(1) ydata(1)]);
                    set(obj.hFig,'WindowButtonMotionFcn',{@obj.MouseAction,'moving_magnifier'})
                    
                case 'moving_magnifier'
                    posM = round(get(obj.hAxesMagnifier,'CurrentPoint'));
                    posM = posM([1 3]) - offset;
                    obj.View(posM(2),posM(1));
                    
                case 'up' % Return to 'mouse_hover' mode
                    set(obj.hFig,'WindowButtonMotionFcn',{@obj.MouseAction,'mouse_hover'})
                    obj.MouseAction([],[],'mouse_hover')
            end
        end
        
        function UpdateFooter(obj,ptrInAxes,posG)
            % update information in the footer
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ptrInAxes
                str = sprintf('%s | Seq %3d | Aln %3d ',obj.seq_names{posG(2)},posG(2),posG(1));
            else
                str = '--- | Seq --- | Aln --- ';
            end
            set(obj.hFooterRight,'String',str)
        end
        
        function UpdatePositions(obj,~,~)
            % update uicontrols positions due to figure size (ResizeFcn)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            posFig = getpixelposition(obj.hMainPanel);
            
            % calculate the maximal possible visible characters
            obj.xWin = max(0,min(obj.x,floor((posFig(3)-obj.SeqLabelExtent-16-10)/obj.cellSize(1))));
            obj.yWin = max(0,min(obj.y,floor((posFig(4)-5-90-20-5-16-5-10-3*obj.cellSize(2)-10)/obj.cellSize(2))));
            
            % check if the main figure is big enough
            if isempty(obj.Data) || obj.xWin<=0 || obj.yWin<=0
                set(obj.hMainPanel,'visible','off')
            else
                set(obj.hMainPanel,'visible','on')
                
                % update settings and positions of both slider
                tmp = obj.x-obj.xWin;
                if tmp > 0
                    set(obj.hSliderHorizontal, ...
                        'max',tmp+1,'value',min(tmp+1,get(obj.hSliderHorizontal,'value')), ...
                        'SliderStep',[1 obj.xWin]/tmp, ...
                        'Position',[1 120 posFig(3) 16])
                else
                    set(obj.hSliderHorizontal, ...
                        'max',1.01,'value',1, ...
                        'sliderstep',[1 inf], ...
                        'Position',[1 120 posFig(3) 16])
                end% if
                
                tmp = obj.y-obj.yWin;
                if tmp>0
                    set(obj.hSliderVertical, ...
                        'max',tmp+1,'value',min(tmp+1,get(obj.hSliderVertical,'value')), ...
                        'SliderStep',[1 obj.yWin]/tmp, ...
                        'Position',[posFig(3)-15 137 16 posFig(4)-137])
                else
                    set(obj.hSliderVertical, ...
                        'max',1.01,'value',1, ...
                        'sliderstep',[1 inf], ...
                        'Position',[posFig(3)-15 137 16 posFig(4)-137])
                end% if
                
                % update axes positions
                set(obj.hAxesConsensus, ...
                    'units','pixel', ...
                    'position',[obj.SeqLabelExtent posFig(4)-obj.cellSize(2)*2-10 [obj.xWin 1].*obj.cellSize], ...
                    'xlim',[0.5 obj.xWin+0.5],'XTick',5:5:obj.xWin, ...
                    'ylim',[0.5 1.5],'yTick',1)
                set(obj.hAxesAlignment, ...
                    'units','pixel', ...
                    'position',[obj.SeqLabelExtent posFig(4)-10-(obj.yWin+3)*obj.cellSize(2) [obj.xWin obj.yWin].*obj.cellSize], ...
                    'xlim',[0.5 obj.xWin+0.5],'XTick',5:5:obj.xWin, ...
                    'ylim',[0.5 obj.yWin+0.5],'yTick',1:1:obj.yWin)
                set(obj.hAxesMagnifier, ...
                    'position',[15 23 posFig(3)-20 100-15])
                
                % update footer
                set(obj.hFooterRight,'position',[4 4 posFig(3)-8 15])
                pos = get(obj.hFooterLeft,'position');
                set(obj.hFooterLeft,'position',[6 6 pos(3) 11])
                
                % prepare the pattern
                obj.seqWin = char(32*ones(2*obj.yWin,2*obj.xWin));
                obj.seqInd = repmat([ true false ; false false ],obj.yWin,obj.xWin);
                
                % update content
                obj.UpdateContent
            end
        end
        
        function UpdateContent(obj,~,~)
            % Update figure content
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~isempty(obj.Data)
                
                % calculate x- and y- index vectors depends on sliders values
                valueH = round(get(obj.hSliderHorizontal,'Value'));
                valueV = round(get(obj.hSliderVertical,'Value'));
                obj.yInd = obj.y-(valueV+obj.yWin-1:-1:valueV)+1;
                % obj.yInd = valueV:valueV+obj.yWin-1;
                obj.xInd = valueH:valueH+obj.xWin-1;
                
                % calculate consensus and score for visible area only!
                % ----------------------------------------------------
                [obj.seqWin(1,obj.seqInd(1,:)),obj.seq_scores] = seqconsensus(obj.seq_txt(:,obj.xInd),'GAPS','ALL');
                obj.seq_scores(obj.seqWin(1,obj.seqInd(1,:))=='-') = max([eps obj.seq_scores(obj.seqWin(1,obj.seqInd(1,:))~='-')]);
                % obj.seqWin(1,obj.seqInd(1,:)) = repmat(' ',1,numel(obj.xInd));
                % obj.seq_scores = ones(1,numel(obj.xInd));
                
                obj.seq_scores = obj.seq_scores/max(obj.seq_scores);
                obj.seq_scores(isnan(obj.seq_scores)) = 0;
                set(obj.hTextConsensus, ...
                    'position',[0 2]+(obj.cellSize./[4 -4]), ...
                    'String',obj.seqWin(1:2,:))
                
                % update consensus score
                % ----------------------
                xx = cell2mat(arrayfun(@(x) {[x x]},1:(numel(obj.seq_scores))+1))-.5;
                yy = [0 1-cell2mat(arrayfun(@(x) {[x x]},obj.seq_scores)) 0]+1.5;
                set(obj.hBarPlotConsensus,'xdata',xx,'ydata',yy)
                
                % update main axes
                obj.seqWin(obj.seqInd) = obj.seq_txt(obj.yInd,obj.xInd);
                set(obj.hTextAlignment, ...
                    'position',[0 2]+(obj.cellSize./[4 -4]), ...
                    'String',obj.seqWin)
                
                % update ticks & labels
                indTick = mod(obj.xInd,5)==0;
                set(obj.hAxesAlignment, ...
                    'Xtick',find(indTick), ...
                    'XTickLabel',obj.xInd(indTick), ...
                    'YTickLabel',obj.seq_names(obj.yInd))
                set(obj.hAxesConsensus, ...
                    'Xtick',find(indTick), ...
                    'XTickLabel',obj.xInd(indTick))
                
                % set magnifier rectangle
                set(obj.hBoxMagnifier, ...
                    'ydata',[obj.yInd(1)-.5 obj.yInd(1)-.5 obj.yInd(end)+.5 obj.yInd(end)+.5], ...
                    'xdata',[obj.xInd(1)-.5 obj.xInd(end)+.5 obj.xInd(end)+.5 obj.xInd(1)-.5])
                
                % update ruler for mgnifier patch
                set(obj.hBoxTopLeftText,'position',[obj.xInd(1)-.5 0 0], ...
                    'String',sprintf('%.0f',obj.xInd(1)-.5))
                set(obj.hBoxTopRightText,'position',[obj.xInd(end)-.5 0 0], ...
                    'String',sprintf('%.0f',obj.xInd(end)-.5))
                set(obj.hBoxLeftTopText,'position',[0 obj.yInd(1)-.5 0], ...
                    'String',sprintf('%.0f',obj.yInd(1)-.5))
                set(obj.hBoxLeftBottomText,'position',[0 obj.yInd(end)-.5 0], ...
                    'String',sprintf('%.0f',obj.yInd(end)-.5))
                
                % consider the edges around the selected area
                % update color background
                switch obj.EnableBackground
                    case 1
                        rgb = get(obj.hBackgroundMagnifier,'CData');
                        set(obj.hBackgroundAlignment,'cdata',rgb(obj.yInd,obj.xInd,:));
                        obj.UpdateSelection
                    case 0
                        obj.Select(obj.xSel,obj.ySel);
                end% switch
                
            end
        end
        
        function UpdateSelection(obj)
            % calculate frame coordinates for the selected area
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~isempty(obj.ySel)
                tmp1 = -floor(obj.SeqLabelExtent/obj.cellSize(1))+1;
                ydata = [min(obj.ySel)-0.5 min(obj.ySel)-0.5 NaN max(obj.ySel)+0.5 max(obj.ySel)+0.5]-obj.yInd(1)+1;
                ydata(find(or(ydata([2 4])<0,ydata([2 4])>obj.yWin+.5))*2)=NaN;
                ydata(ydata>obj.yWin) = obj.yWin+.5;
                ydata(ydata<0) = 0.5;
                set(obj.hSelectRows,'xdata',[tmp1 obj.xWin+0.5 NaN obj.xWin+0.5 tmp1],'ydata',ydata)
            else
                set(obj.hSelectRows,'xdata',NaN,'ydata',NaN)
            end
            if ~isempty(obj.xSel)
                xdata = [min(obj.xSel)-0.5 min(obj.xSel)-0.5 NaN max(obj.xSel)+0.5 max(obj.xSel)+0.5]-obj.xInd(1)+1;
                xdata(find(or(xdata([2 4])<0,xdata([2 4])>obj.xWin+0.5))*2)=NaN;
                xdata(xdata>obj.xWin+0.5) = obj.xWin+.5;
                xdata(xdata<0) = 0.5;
                set(obj.hSelectColumns,'xdata',xdata,'ydata',[-2.5 obj.yWin+0.5 NaN obj.yWin+0.5 -2.5])
            else
                set(obj.hSelectColumns,'xdata',NaN,'ydata',NaN)
            end
            if ~isempty(obj.xSel) && ~isempty(obj.ySel)
                set(obj.hSelectBox,'xdata',xdata([1 1 5 5]),'ydata',ydata([1 5 5 1]))
            else
                set(obj.hSelectBox,'xdata',NaN,'ydata',NaN)
            end
        end
        
    end
end

function FontName = getArchDependentFontName
if ismac
    FontName = 'Monospaced';
elseif ispc
    FontName = 'Lucida Console';
else
    FontName = 'Courier';
end
end
function FontSize = getArchDependentDefaultFontSize
if ismac
    FontSize = 10;
elseif ispc
    FontSize = 8;
else
    FontSize = 10;
end
end
function props = getArchDependentFontProperties(FontName)
% Font extents vary according with architectures and display options on
% every system, getArchDependentFontProperties detects the meaningful
% FontSize for the current font FontName and the width of the edges of a
% text object in the current configuration. 

hf = figure('units','pixels');
ht = text(0,0,repmat('0',1,1),'HorizontalAlign','left',...
    'VerticalAlignment','bottom','FontName',FontName,...
    'FontUnits','pixels','units','pixels');
set(gca,'Units','pixels')
w = zeros(20,10);
h = zeros(20,10);
for fs = 1:20
    for i = 1:10
        set(ht,'string',repmat(randseq(i,'ALPHA','AA'),i,1),'FontSize',fs)
        e = get(ht,'extent');
        w(fs,i) = e(3);
        h(fs,i) = e(4);
    end
end
close(hf)

props.FontSizeReel = find([0;any(diff(w),2)|any(diff(h),2)])';
if ismac
    props.FontSizeReel = props.FontSizeReel(props.FontSizeReel>=5);
end

pw = zeros(20,2);
ph = zeros(20,2);
for fs = 1:20
    pw(fs,:) = polyfit(1:10,w(fs,:),1);
    ph(fs,:) = polyfit(1:10,h(fs,:),1);
end

props.ExtentEdge = [mean(pw(:,2)) mean(ph(:,2))];  
end
