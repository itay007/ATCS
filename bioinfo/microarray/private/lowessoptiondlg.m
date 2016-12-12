function Answer = lowessoptiondlg(DefAns)
% LOWESSOPTIONDLG is derived from INPUTDLG.
% LOWESSOPTIONDLG(OPTIONS), OPTIONS is a cellarray of default setting strings
% OPTIONS{1} - order
% OPTIONS{2} - span
% OPTIONS{3} - robust

% Copyright 2006-2007 The MathWorks, Inc.

Title = 'Options for LOWESS';
WrnTxt = {' (This option can take a long time to calculate.)'};

InputFig=dialog('Visible','off', ...
    'KeyPressFcn'      ,@doFigureKeyPress, ...
    'Name'             ,Title      , ...
    'Pointer'          ,'arrow'    , ...
    'Units'            ,'pixels'   , ...
    'UserData'         ,'Cancel'   , ...
    'Tag'              ,Title      , ...
    'HandleVisibility' ,'callback' , ...
    'NextPlot'         ,'add'      , ...
    'WindowStyle'      ,'normal', ...
    'Resize'           ,'off'       ...
    );

%%% Set Positions %%%
FigWidth=230;
FigHeight=100; %#ok

DefOffset= 5;
BtnWidth  = 53;
BtnHeight = 23;
EditWidth = FigWidth - 2*DefOffset;
EditHeight = 23;

TextInfo.Units              = 'pixels'   ;   
TextInfo.HorizontalAlignment= 'left'     ;
TextInfo.FontSize           = get(0,'FactoryUIControlFontSize');
TextInfo.FontWeight         = get(InputFig,'DefaultTextFontWeight');
TextInfo.HandleVisibility   = 'callback' ;

StInfo=TextInfo;
StInfo.Style              = 'text'  ;

EdInfo=StInfo;
EdInfo = rmfield(EdInfo, 'Style');
EdInfo.BackgroundColor = 'white';
EditStyle = {'popup', 'edit'};

BtnInfo=StInfo;
BtnInfo.Style               = 'pushbutton';
BtnInfo.HorizontalAlignment = 'center';

% Determine # of lines for all Prompts
FigHeight = 6*DefOffset + BtnHeight + 6*EditHeight;
FigPos=get(InputFig,'Position');
FigPos(3)=FigWidth;
FigPos(4)=FigHeight;
set(InputFig, 'Position', FigPos);

YOffset = BtnHeight+3*DefOffset;
EditYOffset = YOffset+DefOffset+2*EditHeight;
QuestYOffset = EditYOffset+(1.5)*EditHeight;

Quest = {'Select order:', 'Enter window size:', 'Use robust fit'};

QuestHandle=[]; %#ok
EditHandle=[];
AxesHandle=axes('Parent',InputFig,'Position',[0 0 1 1],'Visible','off');

%== order string
editStr = {'Linear|Quadratic', DefAns{2}};
for lp=1:2
    QuestHandle(lp)=text('Parent'     ,AxesHandle, ...
        TextInfo     , ...
        'Position'   ,[ DefOffset QuestYOffset+(2-lp)*(DefOffset+2*EditHeight)], ...
        'String'     ,Quest{lp}                 , ...
        'Tag'        ,'Quest' );
    
    EditHandle(lp)=uicontrol(InputFig, ...
        'Style', EditStyle{lp},...
        EdInfo, ...
        'Position'  ,[ DefOffset EditYOffset+(2-lp)*(DefOffset+2*EditHeight) EditWidth EditHeight], ...
        'String'     ,editStr{lp}, ...
        'Tag'        ,'Edit');
end % for lp
set(EditHandle(1), 'Value', str2double(DefAns{1}))

CheckHandle = uicontrol(InputFig,...
    StInfo,...
    'Position', [DefOffset YOffset+0.7*EditHeight EditHeight EditHeight],...
    'Style', 'checkbox',...
    'Tag', 'Check',...
    'Value', DefAns{3});
 
 wrapTextHandle = uicontrol(InputFig, ...
     TextInfo, ...
     'Position',[DefOffset+20 YOffset EditWidth-23 1.5*EditHeight], ...
     'Style', 'text',...
     'Tag','Check'); %#ok
 
 wrapTxt = textwrap(wrapTextHandle, {strcat(Quest{3}, WrnTxt{:})});
 set(wrapTextHandle, 'String', wrapTxt);
              
OKHandle=uicontrol(InputFig     ,              ...
                   BtnInfo      , ...
                   'Position'   ,[ FigWidth-2*BtnWidth-2*DefOffset DefOffset BtnWidth BtnHeight ] , ...
                   'KeyPressFcn',@doControlKeyPress , ...
                   'String'     ,'OK'        , ...
                   'Callback'   ,@doCallback , ...
                   'Tag'        ,'OK'        , ...
                   'UserData'   ,'OK'          ...
                   );
               
CancelHandle=uicontrol(InputFig, ...
                       BtnInfo      , ...
                       'Position'   ,[ FigWidth-BtnWidth-DefOffset DefOffset BtnWidth BtnHeight ]           , ...
                       'KeyPressFcn',@doControlKeyPress            , ...
                       'String'     ,'Cancel'    , ...
                       'Callback'   ,@doCallback , ...
                       'Tag'        ,'Cancel'    , ...
                       'UserData'   ,'Cancel'      ...
                       ); %#ok
                   
setDefaultButton(InputFig, OKHandle)
set(InputFig,'Visible','on');
drawnow;
uiwait(InputFig);

if ishandle(InputFig)
    Answer={};
    if strcmp(get(InputFig,'UserData'),'OK'),
        Answer=cell(3,1);
        Answer{1} = get(EditHandle(1), 'Value');
        Answer(2)=get(EditHandle(2),{'String'});
        Answer{3} = get(CheckHandle, 'Value');
    end
    delete(InputFig);
else
    Answer={};
end

%------------------------------------------------------------%
function doFigureKeyPress(obj, evd) %#ok
switch(evd.Key)
 case {'return','space'}
  set(gcbf,'UserData','OK');
  uiresume(gcbf);
 case {'escape'}
  delete(gcbf);
end

function doControlKeyPress(obj, evd) %#ok
switch(evd.Key)
 case {'return'}
  if ~strcmp(get(obj,'UserData'),'Cancel')
      set(gcbf,'UserData','OK');
      uiresume(gcbf);
  else
      delete(gcbf)
  end
 case 'escape'
  delete(gcbf)
end

function doCallback(obj, evd) %#ok
if ~strcmp(get(obj,'UserData'),'Cancel')
    set(gcbf,'UserData','OK');
    uiresume(gcbf);
else
    delete(gcbf)
end

function setDefaultButton(figHandle, btnHandle) %#ok
% First get the position of the button.
btnPos = getpixelposition(btnHandle);

% Next calculate offsets.
leftOffset   = btnPos(1) - 1;
bottomOffset = btnPos(2) - 2;
widthOffset  = btnPos(3) + 3;
heightOffset = btnPos(4) + 3;

% Create the default button look with a uipanel.
% Use black border color even on Mac or Windows-XP (XP scheme) since
% this is in natve figures which uses the Win2K style buttons on Windows
% and Motif buttons on the Mac.
h1 = uipanel(get(btnHandle, 'Parent'), 'HighlightColor', 'black', ...
    'BorderType', 'etchedout', 'units', 'pixels', ...
    'Position', [leftOffset bottomOffset widthOffset heightOffset]);

% Make sure it is stacked on the bottom.
uistack(h1, 'bottom');

                 