function bioinfoFeedback(componentName)
%BIOINFOFEEDBACK Open a mail page with for feedback about the component.

tlbx = ver('bioinfo');

mailstr = ['mailto:bioinfo-feedback@mathworks.com?subject=',...
    'Feedback%20for%20' componentName '%20in%20Bioinformatics',...
    '%20Toolbox%20',tlbx(1).Version];
web(mailstr)
end