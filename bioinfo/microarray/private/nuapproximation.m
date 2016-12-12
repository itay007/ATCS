function NU = nuapproximation
% NUAPPROXIMATION numerical approximation of the special function nu.

% HYBRIDNU = NUAPPROXIMATION returns an anonymous function NU for
% numerical approximation for the special function nu in compuattion for
% the tail probabilities of the statistic T_2+. (The T_2+ is the same as
% the statistic Z_3 in the reference paper.)
% v(x) = 2(x^-2)exp{-2sigma(l^-1 phy(-1/2*(xl)^1/2))}
% 
% References:
% [1] Q. Yao, 'Tests for Change-Points with Epidemic Alternatives',
%     Biometrika (1993), 80, 179-191.  

%   Copyright 2007-2009 The MathWorks, Inc. 


% Define the summary terms
ft = @(m,x) 1./m .*normcdf(-sqrt(m).*x/2);

% Set the tolerance in numerical approximation
errortarget = 1e-8;

% Approximate the number of terms in summary for different x
xmin = 0.001; % minimum value x can take
xmax = 20;    % maximum value x can take
xnum = 10000; % resolution for search space in x

% Search space in x is log to enhance resolution
x = exp(log(xmin):(log(xmax)-log(xmin))/(xnum-1):log(xmax));
y = zeros(numel(x),1);

for i = 1:numel(x)
    le = 1;
    ri = 1e9;
    step = round((ri-le)/10);     
    fta = ft(1,x(i))*errortarget;
    while step
        q = find(ft(le+(0:10)*step,x(i))<=fta,1);
        ri = le + step*(q-1);
        le = le + step*(q-2);
        step = round((ri-le)/10);
    end
    y(i) = le;    
end

h = [true;diff(y)~=0];
x = x(h);
y = y(h);

% Create an anonymous function for computing number of terms 
nterms = @(z) y(find(x>z,1));

% It is sufficient to use approximation for small x (see ref.) Everything
% below 0.003 seems to have round-off error, you can change xmin to .001 or
% smaller if you want to see it.
xmin = 0.0046; % minimum value x can take 
xmax = 1000000;    % maximum value x can take
xnum = 999; % resolution for search space in x
% In logarithm space
x = [0 exp(log(xmin):(log(xmax)-log(xmin))/(xnum-1):log(xmax))]'; 

v = zeros(numel(x),1);
for i=1:numel(x)
  v(i) = 2/x(i)/x(i)*exp(-2*sum(ft(1:max(10,nterms(x(i))),x(i))));
  %disp(sprintf('%d %d %f %f',i,nterms(x(i)),x(i),v(i)))
end

% Find the slope for the line that departs from (0,1) and touches v just at
% one sigle point:
[mm,mmh] = max((1-v)./x);
vc = v;
vc(1:mmh) = 1-x(1:mmh).*mm;  % 0.575240945666038 (0.583 in Ref.)

% Create an anonymous function for computing v on the fly...
NU = @(X) interp1q(x,vc,X);

