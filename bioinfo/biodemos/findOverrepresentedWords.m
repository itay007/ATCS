function [nmersSorted, freqDiffSorted] = findOverrepresentedWords(seq, seq0, W)
% FINDOVERREPRESENTEDWORDS helper for evemotifdemo

% Copyright 2007 The MathWorks, Inc.

%=== find and count words of length W
nmers0 = nmercount(seq0, W);
nmers = nmercount(seq, W);

%=== compute frequency of words 
f = [nmers{:,2}]/(length(seq) - W + 1);
f0 = [nmers0{:,2}]/(length(seq0) - W + 1);

%=== determine words common to both set 
[nmersInt, i1, i2] = intersect(nmers(:,1),nmers0(:,1));
freqDiffInt = (f(i1) - f0(i2))';

%=== determine words specific to one set only
[nmersXOr, i3, i4] = setxor(nmers(:,1),nmers0(:,1));
c0 = nmers(i3,1);
d0 = nmers0(i4,1);
nmersXOr = [c0; d0]; 
freqDiffXOr = [f(i3) -f0(i4)]';

%=== define all words and their difference in frequency (margin)
nmersAll = [nmersInt; nmersXOr];
freqDiff = [freqDiffInt; freqDiffXOr];

%=== sort according to descending difference in frequency
[freqDiffSorted, freqDiffSortedIndex] = sort(freqDiff, 'descend'); 
nmersSorted = nmersAll(freqDiffSortedIndex);

