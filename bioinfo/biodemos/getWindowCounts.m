function [counts,pval,fdr] = getWindowCounts(bm,t,w,l)
%GETWINDOWCOUNTS finds windows with significant coverage.
%
% GETWINDOWCOUNTS is a helper function for MBDSEQDEMO.
%
% [COUNTS,PVAL,FDR] = GETWINDOWCOUNTS(BM,T,W,L) returns the counts for
% each window and fits a right-truncated negative binomial as a null
% distribution to find windows where the counts are unlikely to be by
% chance. T is the truncating threshold. W is a vector indicating the
% start of each window. L is a scalar indicating the length of the
% windows.

% Copyright 2012 The MathWorks, Inc.

rtnbinpdf = @(x,p1,p2,t) nbinpdf(x,p1,p2) ./ nbincdf(t-1,p1,p2);
rtnbinfit = @(x,t) mle(x,'pdf',@(x,p1,p2) rtnbinpdf(x,p1,p2,t),'start',nbinfit(x),'lower',[0 0]);

counts = getCounts(bm,w,w+l-1,'independent',true,'overlap','start');
pn = rtnbinfit(counts(counts<t),t);
pval = 1 - nbincdf(counts,pn(1),pn(2));
fdr = mafdr(pval,'bhfdr',true);
