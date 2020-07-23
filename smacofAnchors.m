% Copyright (c) 2020
% Author: Carmelo Di Franco 
% Email: carmelodfr@gmail.com

% A portion of this code is copyrighted by:

% TOSCA = Toolbox for Surface Comparison and Analysis
% Web: http://tosca.cs.technion.ac.il
% Version: 1.0
%
% (C) Copyright Michael Bronstein, 2008
% All rights reserved.
%
% License:
%
% ANY ACADEMIC USE OF THIS CODE MUST CITE THE ABOVE REFERENCES. 
% ANY COMMERCIAL USE PROHIBITED. PLEASE CONTACT THE AUTHORS FOR 
% LICENSING TERMS. PROTECTED BY INTERNATIONAL INTELLECTUAL PROPERTY 
% LAWS. PATENTS PENDING.
function [X,hist] = smacofAnchors(D,X0,W,n_anchors,MDSparam,fixedheight)

% check input correctness
if nargin < 5,
    error('Incorrect number of arguments, exiting.')     
end

iter     = MDSparam.iter;
verbose  = MDSparam.verbose;
xhistory = MDSparam.xhistory;
rtol     = MDSparam.rtol;
atol     = MDSparam.atol;

if size(D,1) ~= size(D,2),
    error('Matrix D must be square, exiting.')     
end

if size(D,1) ~= size(X0,1),
    error('X0 and D dimensions mismatch, exiting. X0 must be a size(D,1)*dim matrix.')     
end

if rtol < 0
    error('rtol must be non-negative, exiting.')     
end


% check input and output flags
if strcmp(lower(verbose),'iter'),
    VERBOSE = 1;
else
    VERBOSE = 0;
end

if nargout == 2,
    HISTORY = 1;
else
    HISTORY = 0;
end

if strcmp(lower(xhistory),'on'),
    XHISTORY = 1;
else
    XHISTORY = 0;
end


% initialize
iii = 1;
Z = X0;
X = X0;
D_ = calc_D (X);

% initialize history
if HISTORY | VERBOSE,
    hist.time = zeros(iter);
    %hist.s = zeros(iter);
    hist.s(1) = calc_stress(X0,D,W);
end

if XHISTORY
   hist.X{1} = X0; 
end

if HISTORY & VERBOSE,
    fprintf(1,'iter         stress   time (sec)\n') 
    fprintf(1,'INIT   %12.3g   ----------\n', hist.s(1)) 
end


nodes = size(D,1);
n = size(D,1) - n_anchors;

V = calc_V(W);
V11 = V(1:n,1:n);
V12 = V(1:n,n+1:nodes);
%V21 = V(n+1:nodes,1:n);
%V22 = V(n+1:nodes,n+1:nodes);

Vp11 = pinv(V11);
    
    
    
while (iii <= iter),
    t = cputime;       

     
    B_ = calc_B(D_,D,W);
    B11 = B_(1:n,1:n);
    B12 = B_(1:n,n+1:nodes);
   % B21 = B_(n+1:nodes,1:n);
   % B22 = B_(n+1:nodes,n+1:nodes);
    
    Z1 = Z(1:n,:);
    Z2 = Z(n+1:nodes,:);
    X2 = X0(n+1:nodes,:);

 
     
    X1 = Vp11*( + B11*Z1 +  B12*Z2 - V12*X2);
    X = [ X1; X2];
 
    if (fixedheight == 1 && size(X0,2) == 3)
        X(:,3) = X0(:,3);
        % it seems an override, but it is exactly the same as using only 2
        % columns (x,y) and do not use z since they are not affected. see
        % journal paper to see the math.
        % it should be X1(1:2,:)  = Vp11*( + B11*Z1(1:2,:)  +  B12*Z2(1:2,:)  - V12*X2(1:2,:) );
    end
    
    D_ = calc_D (X);
    S = calc_S (D,D_,W);
    Z = X;
    
    %hist.X{iii} = X; 
    %hist.time(iii) = cputime-t;
    hist.s(iii) = calc_stress(X,D,W);
     
    
    if VERBOSE
        fprintf(1,'%4d   %12.3g   %10.3g\n', iii,hist.s(iii),hist.time(iii)) 
    end
    
    
    % check stopping conditions
    if S < atol
        fprintf(1,'atol=%g reached, exiting, iters = %g,\n',atol,iii)
        return 
    end

    if (iii > 1)
        %if (hist.s(iii-1)/hist.s(iii)-1) < rtol,
        if ( abs(hist.s(iii-1) - hist.s(iii))) < rtol
        
            fprintf(1,'rtol=%g reached, exiting, iters = %g,\n',rtol,iii);
            return 
        end
    end

    iii = iii+1;
        
end    



% SERVICE FUNCTIONS

% compute the stress 
function [S] = calc_stress (X,D,W)
D_ = calc_D (X);
S  = calc_S (D,D_,W);
return

function [S] = calc_S (D,D_,W)
d = triu(W.*(D - D_).^2,1);
S = sum(d(:));
return

function [D] = calc_D (X)
D = zeros(size(X,1));
for k=1:size(X,1),
    xk = repmat(X(k,:),size(X,1),1);
    D(:,k) = sqrt(sum((X - xk).^2, 2));
end    
return;

function [B] = calc_B (D_,D,W)
 

B = zeros(size(D));
i = find(D_(:) ~= 0);
B(i) = - (W(i).*D(i))./D_(i);
B = B - diag(diag(B));
d = sum(B);
B = B - diag(d);
 
return;


function [V] = calc_V (W)
 
nodes = size(W,1);

V = zeros(size(W));
i = find(W(:) ~= 0);
V(i) = - W(i);
V = V - diag(diag(V));
d = sum(V);
V = V - diag(d);
return;
