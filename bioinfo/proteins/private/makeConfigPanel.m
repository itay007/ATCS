function [panel,windowcombobox,edgeweightcombospinner,weightgroup,buttons] = makeConfigPanel(handles)
% MAKECONFIGPANEL Utility function for PROTEINPLOT
%
% function [panel,windowcombobox,edgeweighttextfield] = makeConfigPanel(handles)
% create configuration panel

%   Copyright 2003-2009 The MathWorks, Inc.


import javax.swing.*; import java.awt.*; import java.awt.event.*;
import com.mathworks.mwswing.*;

% create panel
panel = MJPanel([]);
% panel = MJPanel(GridLayout(3,2));

x = 25;
y = 20;

w = 100;
h = 20;

windowlabel = MJLabel('Window Size:', [], JLabel.LEFT);
windowcombobox = MJComboBox(num2cell(handles.windowrange)); 
h_windowcombobox = handle(windowcombobox,'callbackproperties');
edgeweightlabel = MJLabel('Edge Weight:', [], JLabel.LEFT);

spinnerModel = javax.swing.SpinnerNumberModel(1.00, 0, 1.00, 0.01);
edgeweightcombospinner = com.mathworks.mwswing.MJSpinner(spinnerModel); 

% % edgeweightcombospinner = com.mathworks.page.utils.ComboSpinner(1,0,1,.01); 
h_edgeweightcombospinner = handle(edgeweightcombospinner,'callbackproperties');


windowlabel.setBounds(x,y,w,h);
panel.add(windowlabel);


windowcombobox.setBounds(3*x + w,h,w,h);
% windowcombobox.setPreferredSize(Dimension(75,25));
% windowcombobox.setSize(Dimension(75,25));
windowcombobox.setSelectedItem(windowcombobox.getItemAt(find(handles.windowrange == handles.windowlength) - 1));
set(h_windowcombobox,'ActionPerformedCallback', {@windowsizechanged,windowcombobox,handles.proteinplotfig});
panel.add(windowcombobox);


edgeweightlabel.setBounds(x,2*y + h,w,h);
edgeweightlabel.setToolTipText('Values should be between 0 and 1');
panel.add(edgeweightlabel);


edgeweightcombospinner.setBounds(3*x + w,2*y + h,w,h);
% edgeweightcombospinner.setSize(Dimension(75,25));
% edgeweightcombospinner.setInsets(Insets(20,0,20,0));

% % set(h_edgeweightcombospinner,'ActionPerformedCallback', {@edgeweightchanged,edgeweightcombospinner,handles.proteinplotfig});
set(h_edgeweightcombospinner,'StateChangedCallback', {@edgeweightchanged,edgeweightcombospinner,handles.proteinplotfig});
set(h_edgeweightcombospinner,'FocusLostCallback', {@edgeweightchanged,edgeweightcombospinner,handles.proteinplotfig});
panel.add(edgeweightcombospinner);


smoothinglabel = MJLabel('Smoothing:');
smoothinglabel.setBounds(x,3*y + 2*h,w,h);
panel.add(smoothinglabel);

buttonpanel = MJPanel(GridLayout(3,1));

linearradio = MJRadioButton('Linear'); h_linearradio = handle(linearradio,'callbackproperties');
expradio = MJRadioButton('Exponential'); h_expradio = handle(expradio,'callbackproperties');
lowessradio = MJRadioButton('Lowess'); h_lowessradio = handle(lowessradio,'callbackproperties');

buttons.linear = linearradio;
buttons.exp = expradio;
buttons.lowess = lowessradio;



linearradio.setSelected(true);
% linearradio.setBounds(25,25,buttonwidth,25);
linearradio.setSize(w,h);
set(h_linearradio,'ActionPerformedCallback', {@linearradioselected,linearradio,handles.proteinplotfig,edgeweightcombospinner});
buttonpanel.add(linearradio);


% expradio.setBounds(25,50,buttonwidth,25);
expradio.setSize(w,h);
set(h_expradio,'ActionPerformedCallback', {@expradioselected,expradio, handles.proteinplotfig,edgeweightcombospinner});
buttonpanel.add(expradio);


% lowessradio.setBounds(25,75,buttonwidth,25);
lowessradio.setSize(w,h);
set(h_lowessradio,'ActionPerformedCallback', {@lowessradioselected,lowessradio,handles.proteinplotfig,edgeweightcombospinner});
buttonpanel.add(lowessradio);

weightgroup = ButtonGroup;
weightgroup.add(linearradio);
weightgroup.add(expradio);
weightgroup.add(lowessradio);

buttonpanel.setBounds(3*x + w,3*y + 2*h,w,3*h );

buttonpanel.setBorder(...
    BorderFactory.createCompoundBorder(...
                                          BorderFactory.createLineBorder(Color.gray),...
                                          BorderFactory.createEmptyBorder(5,5,5,5)));

panel.add(buttonpanel);
% panel.setSize(275,350);



function windowsizechanged(jcb,data,h,hfig) %#ok
handles = guidata(hfig);
handles.windowlength = handles.windowrange(h.getSelectedIndex + 1);
guidata(handles.proteinplotfig,handles);

if handles.autoapply
    proteinplot('analyze', h, [], handles)
end

function edgeweightchanged(jtf,data,h,hfig) %#ok
handles = guidata(hfig);
% % handles.edgeweight = str2num(h.getText); %#ok
handles.edgeweight = h.getValue; %#ok
guidata(handles.proteinplotfig,handles);
if handles.autoapply
    proteinplot('analyze', h, [], handles)
end

function linearradioselected(jb,data,h,hfig,ewcs) %#ok
handles = guidata(hfig);
handles.uselinear = true;
handles.useexp = false;
handles.uselowess = false;
guidata(handles.proteinplotfig,handles);
awtinvoke(ewcs, 'setEnabled(Z)', true);
if handles.autoapply
    proteinplot('analyze', h, [], handles)
end

function expradioselected(jb,data,h,hfig,ewcs) %#ok
handles = guidata(hfig);
handles.useexp = true;
handles.uselinear = false;
handles.uselowess = false;
guidata(handles.proteinplotfig,handles);
awtinvoke(ewcs, 'setEnabled(Z)', true);
if handles.autoapply
    proteinplot('analyze',h, [], handles)
end

function lowessradioselected(jb,data,h,hfig,ewcs) %#ok
handles = guidata(hfig);
handles.uselowess = true;
handles.uselinear = false;
handles.useexp = false;
guidata(handles.proteinplotfig,handles);
awtinvoke(ewcs, 'setEnabled(Z)', false);
if handles.autoapply
    proteinplot('analyze', h, [], handles)
end
