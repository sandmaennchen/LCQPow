function [z, x, b, y, stats] = miqpsol(miqp,params)
%MIQPSOL solver for MIQPs using LCQPow.
%
%
%                min   1/2*[z; x]'H[z; x] + [z; x]'g
%                s.t.     0 <=  z <= ubz
%                       lbx <=  x <= ubx      {optional}
%                       lbA <= Ax <= ubA      {optional}
%                               z integer
%
% The struct miqp is expected to provide the following fields:
% Mandatory fields
%   nz, nx: number of integer, continuous optimization variables
%   H:
%   ubz:
%
% Optional fields:
%   g:
%   lbx, ubx: lower and upper bounds on x
%   A, lbA, ubA:

mandatory_fields = {'nz', 'nx', 'H', 'ubz'};
optional_fields = {'g', 'lbx', 'ubx', 'A', 'lbA', 'ubA'};

% check if mandatory fields are present
for k=1:length(mandatory_fields)
    f = mandatory_fields{k};
    if ~isfield(miqp, f)
        error(['Mandatory field ' f ' is missing!']);
    end
end

nx = miqp.nx;
nz = miqp.nz;

ubz = miqp.ubz;
H = miqp.H;

if ~(isvector(ubz) && size(ubz, 1) == nz)
    error(['ubz has wrong size, expected vector of dimension nz = ' num2str(nz)]);
end

if ~(ismatrix(H) && size(H, 1) == nx+nz && size(H, 2) == nx+nz)
    error(['H has wrong size, expected quadratic matrix of dimension nz + nx = ' num2str(nx+nz)]);
end

% check which optional fields are present
for k=1:length(optional_fields)
    f = optional_fields{k};
    
    switch f
        case 'g'
            if isfield(miqp, 'g')
                g = miqp.g;
                if ~(isvector(g) && size(g, 1) == nz+nx)
                    error(['g has wrong size, expected vector of dimension nz + nx = ' num2str(nx+nz)]);
                end
            else
                g = zeros(nx+nz, 1);
            end
        case 'lbx'
            if isfield(miqp, 'lbx')
                lbx = miqp.lbx;
                if ~(isvector(lbx) && size(lbx, 1) == nx)
                    error(['lbx has wrong size, expected vector of dimension nx = ' num2str(nx)]);
                end
            else
                lbx = -inf*ones(nx, 1);
            end
        case 'ubx'
            if isfield(miqp, 'ubx')
                ubx = miqp.ubx;
                if ~(isvector(ubx) && size(ubx, 1) == nx)
                    error(['ubx has wrong size, expected vector of dimension nx = ' num2str(nx)]);
                end
            else
                ubx = inf*ones(nx, 1);
            end
        case 'A'
            if isfield(miqp, 'A')
                A = miqp.A;
                if ~(ismatrix(A) && size(A, 2) == nz+nx)
                    error(['A has wrong size, expected matrix of dimension nz + nx = ' num2str(nz+nx)]);
                end
                n_ineq = size(A, 1);
            else
                n_ineq = 0;
            end
        case 'lbA'
            if isfield(miqp, 'lbA')
                lbA = miqp.lbA;
                if ~(isvector(lbA) && size(lbA, 1) == n_ineq)
                    error(['lbA has wrong size, expected vector of dimension n_ineq = ' num2str(n_ineq)]);
                end
            else
                if n_ineq > 0
                    lbA = -inf*ones(n_ineq, 1);
                end
            end
        case 'ubA'
            if isfield(miqp, 'ubA')
                ubA = miqp.ubA;
                if ~(isvector(ubA) && size(ubA, 1) == n_ineq)
                    error(['ubaA has wrong size, expected vector of dimension n_ineq = ' num2str(n_ineq)]);
                end
            else
                if n_ineq > 0
                    ubA = inf*ones(n_ineq, 1);
                end
            end
    end
end

% number of binary variables per integer variable
n_comp = ceil(log2(ubz + 1));

% total number of binary variables
nb = sum(n_comp);

% compute powers of two needed to represent integer variables
max_n_comp = max(n_comp);
powers = ones(max_n_comp, 1);
for i=1:max_n_comp+1
    powers(i+1) = 2*powers(i);
end

% build equality constraints matrix that maps
% binary variables to integer variables
C = zeros(nz, nb);
offset = 1;

for i=1:nz
    C(i, offset:(offset+n_comp(i)-1)) = powers(1:n_comp(i));
    offset = offset + n_comp(i);
end

C_tilde = [C, -eye(nz), zeros(nz, nx), zeros(nz, 1)];

disp("number of binary variables per integer variable");
disp(n_comp);
disp("C");
disp(C_tilde);

%%% build LCQP %%%

% number of optimization variables
nw = nb + nz + nx + 1;

% cost
H_tilde = zeros(nw, nw);
H_tilde(nb+1:nb+nz+nx, nb+1:nb+nz+nx) = H;

g_tilde = zeros(nw, 1);
g_tilde(nb+1:nb+nz+nx) = g;

% constraints
if n_ineq > 0
    A = [zeros(n_ineq, nb), A, zeros(n_ineq, 1)];
    C_tilde = [C_tilde; A];
end

lbC_tilde = zeros(nz+n_ineq, 1);
ubC_tilde = zeros(nz+n_ineq, 1);

if n_ineq > 0
    lbC_tilde(nz+1:nz+n_ineq) = lbA;
    ubC_tilde(nz+1:nz+n_ineq) = ubA;
end

% lower bound (TODO: 0 for binary and/or integer variables?)
lb = -inf*ones(nw, 1);
lb(nb+nz+1:nb+nz+nx) = lbx;
lb(end) = 1;

% upper bound (TODO: 1 for binary variables?)
ub = inf*ones(nw, 1);
ub(nb+1:nb+nz) = ubz;
ub(nb+nz+1:nb+nz+nx) = ubx;
ub(end) = 1;

S1 = [eye(nb), zeros(nb, nz+nx+1)];
S2 = [-eye(nb), zeros(nb, nz+nx), ones(nb, 1)];

params.x0 = [zeros(nb, 1); params.z0; params.x0; ones(1,1)];

[w, y, stats] = LCQPow(H_tilde,g_tilde,S1,S2,C_tilde,lbC_tilde,ubC_tilde,lb,ub,params);

b = w(1:nb);
z = w(nb+1:nb+nz);
x = w(nb+nz+1:nb+nz+nx);

disp("z");
disp(z);

disp("x");
disp(x);

disp("b");
disp(b);
end

