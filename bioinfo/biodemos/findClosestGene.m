function [s,msg] = findClosestGene(a,w)
%FINDCLOSESTGENE finds genes and promoter regions that overlap with range.
%
% FINDCLOSESTGENE is a helper function for MBDSEQDEMO.
%
% [S,MSG] = FINDCLOSESTGENE(A,W) returns a subset of the GFFAnnotation
% object with the gene that overlaps the region specified by the two
% element vector W. MSG is a string indicating if the region is intergenic,
% intragenic, or in a proximal promoter region (defined as -500/100 5'
% region from the gene's TSS).

%  Copyright 2012 The MathWorks, Inc.

d = getData(a);
geneStart = [d.Start];
geneStop = [d.Stop];
geneDir = [d.Strand]=='+';

promoterStart(geneDir) = geneStart(geneDir) - 500;
promoterStop(geneDir) = geneStart(geneDir) + 100;
promoterStart(~geneDir) = geneStop(~geneDir) - 100;
promoterStop(~geneDir) = geneStop(~geneDir) + 500;

r = find((w(2)>=promoterStart) & (w(1)<=promoterStop));
if ~isempty(r)
    msg = sprintf('Prox. Promoter (%s)',d(r(1)).Feature);
    s = GFFAnnotation(d(r(1)));
else
    r = find((w(2)>=geneStart) & (w(1)<=geneStop));
    if ~isempty(r)
        msg = sprintf('Intergenic (%s)',d(r(1)).Feature);
        s = GFFAnnotation(d(r(1)));
    else
        msg = 'Intragenic';
        s = GFFAnnotation(d([]));
    end
end
