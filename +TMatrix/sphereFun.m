function varargout = sphereFun(kind,n,x,m)
%% sphereFun - radial and angular partial-wave functions, with derivatives
%
%   [Z0,Z1,Z2,Z3] = sphereFun('j',n,x)    spherical Bessel j_n and d/dx up to 3rd order
%   [Z0,Z1,Z2,Z3] = sphereFun('h',n,x)    spherical Hankel h_n^(1), same
%   [L0,L1,L2,L3] = sphereFun('P',n,u)    Legendre P_n and d/du up to 3rd order
%
% Radial derivatives are taken with respect to the ARGUMENT x = k*r; the
% caller multiplies by k^m to get d^m/dr^m.

switch kind
    case {'j','h'}
        % z_n(x) from besselj/bessely of half-integer order
        z = @(nn) sph(kind,nn,x);
        % derivative: z_n'(x) = z_{n-1}(x) - (n+1)/x * z_n(x), applied repeatedly
        % (more stable than differentiating the three-term recurrence)
        Z = zeros(4,numel(x));
        for m0 = 0:3
            Z(m0+1,:) = dz(kind,n,x,m0,z);
        end
        varargout = {Z(1,:),Z(2,:),Z(3,:),Z(4,:)};

    case 'P'
        % P_n and derivatives from the recurrence differentiated m times:
        %   (k+1) P_{k+1}^(m) = (2k+1)[ u P_k^(m) + m P_k^(m-1) ] - k P_{k-1}^(m)
        u = x(:).';
        nd = 5;                       % orders 0..4 (the 4th is needed for P_n^1)
        Pm = zeros(nd,n+2,numel(u));  % Pm(m+1,k+1,:) = P_k^(m)
        Pm(1,1,:) = 1;
        if n >= 0
            Pm(1,2,:) = u;
            Pm(2,2,:) = 1;
        end
        for k = 1:n
            for m0 = 0:nd-1
                term = (2*k+1)*(u.*squeeze(Pm(m0+1,k+1,:)).');
                if m0 >= 1
                    term = term + (2*k+1)*m0*squeeze(Pm(m0,k+1,:)).';
                end
                Pm(m0+1,k+2,:) = (term - k*squeeze(Pm(m0+1,k,:)).')/(k+1);
            end
        end
        varargout = {squeeze(Pm(1,n+1,:)).', squeeze(Pm(2,n+1,:)).', ...
                     squeeze(Pm(3,n+1,:)).', squeeze(Pm(4,n+1,:)).', ...
                     squeeze(Pm(5,n+1,:)).'};
    otherwise
        error('kind must be j, h or P');
end
end

function v = sph(kind,n,x)
% spherical Bessel/Hankel from the cylindrical ones of half-integer order
s = sqrt(pi./(2*x));
switch kind
    case 'j', v = s.*besselj(n+0.5,x);
    case 'h', v = s.*(besselj(n+0.5,x) + 1i*bessely(n+0.5,x));
end
end

function v = dz(kind,n,x,m,~)
% m-th derivative of z_n(x) with respect to x, by repeated application of
%   z_n'(x) = z_{n-1}(x) - (n+1)/x z_n(x)
% expanded up to m=3 in terms of z_{n-3..n}
if m == 0
    v = sph(kind,n,x); return
end
% first derivative of an arbitrary combination: the result is represented as
% a sum of coefficients times z_{n-k}. Starts with z_n.
coef = {@(x) ones(size(x))};   % coefficient of z_{n-0}
order = 0;
for it = 1:m
    newCoef = cell(1,numel(coef)+1);
    [newCoef{:}] = deal(@(x) zeros(size(x)));
    for k = 1:numel(coef)
        ck = coef{k}; nk = n - order(k);
        % d/dx [ c * z_{nk} ] = c' z_{nk} + c*( z_{nk-1} - (nk+1)/x z_{nk} )
        dc = @(x) derivNum(ck,x);
        newCoef{k}   = sumFun(newCoef{k}, @(x) dc(x) - (nk+1)./x.*ck(x));
        newCoef{k+1} = sumFun(newCoef{k+1}, ck);
    end
    coef = newCoef;
    order = 0:numel(coef)-1;
end
v = zeros(size(x));
for k = 1:numel(coef)
    v = v + coef{k}(x).*sph(kind,n-(k-1),x);
end
end

function f = sumFun(a,b), f = @(x) a(x) + b(x); end

function d = derivNum(f,x)
% derivative of the rational coefficients in x (powers of 1/x): high-order
% central difference, accurate enough for these smooth, cheap functions
h = 1e-5*max(abs(x),1);
d = (-f(x+2*h) + 8*f(x+h) - 8*f(x-h) + f(x-2*h))./(12*h);
end
