function [m0,m1,s0,s1,p] = em2c1d(x)
%EM2C1D provides an expectation-maximization estimation
%
% [m0,m1,s0,s1,p] = EM2C1D(x) utilizes an expectation-maximization
% algorithm to estimate a mixture of two Gaussians in one dimension. 

%   Copyright 2003-2010 The MathWorks, Inc.


% References: "Signal Background Estimation and Baseline Correction
% Algorithms for Accurate DNA Sequencing", L. Andrade, E. Manolakos

if numel(x)==1
    m0 = x;
    m1 = x;
    s0 = 0;
    s1 = 0;
    p = 0.5;
    return
end

% Initial safe guesses
m0=quantile(x,.10);
m1=quantile(x,.90);
s0=std(x); s1=s0;
p=0.5;

% Change in log-likelihood needed to call it a local optimum and leave the
% loop
changeInLThershold = .01; % i.e. 1% 

L = -inf;
stayInLoop = true;

while stayInLoop
    [cid,newL] = expectation(x,m0,m1,s0,s1,p);
    [m0,m1,s0,s1,p] = maximization(x,cid);

    stayInLoop = (newL > L) && ((newL - L)/newL < changeInLThershold );
    L = newL;
end

function [cid,L] = expectation(x,m0,m1,s0,s1,p)
p0 = normpdf(x,m0,s0) * p;
p1 = normpdf(x,m1,s1) * (1-p);
sp = p0 + p1;
cid = p0 < p1;

% Make sure that the two classes are represented with at least one sample
% even if the posterior < 0.5
if all(cid)
   [~,h]=min(log(p1)-log(sp));
   cid(h)=0;
elseif all(~cid)
   [~,h]=min(log(p0)-log(sp));
   cid(h)=1;
end
    
% log-likelihood of the model
L = sum(log(p0(~cid))) + sum(log(p1(cid))) - sum(log(sp));

function [m0,m1,s0,s1,p] = maximization(x,cid)
m0 = mean(x(~cid));
m1 = mean(x(cid));
s0 = std(x(~cid));
s1 = std(x(cid));
p = sum(~cid)/numel(x);
