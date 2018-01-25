% This will be a system which is homogeneous, steady state. It will be 
% modelled in three-dimensions as a test. There will be a pre-determined
% pressure coming into one side of the system and coming out the other
% side such that Darcy flow is achieved. A cross section of this will be
% compared to the analytical solution for Darcy flow. Flow is in the z
% direction.
%
% Rows (R) go in the x-direction, columns (C) in the y-direction, and
% layers (L) in the z-direction. Indices start at 1, end at R, C, and L.
% The model's indices are numbered from left to right (x), top to bottom
% (y), front to back (z), such that, for node i in x, j in y, and k in z,
% the center node is R*C*(k-1)+C*(j-1)+i. 
%
% The indices which are used are labeled below
%
%    IC = R*C*(k-1)+C*(j-1)+i;   IC is the index for the center node
%    PZ = IC + R*C;              PZ is the index for +1 in the Z direction
%    PY = IC + C;                PY is the index for +1 in the Y direction
%    PX = IC + 1;                PX is the index for +1 in the X direction
%    MZ = IC - R*C;              MZ is the index for -1 in the Z direction
%    MY = IC - C;                MY is the index for -1 in the Y direction
%    MX = IC - 1;                MX is the index for -1 in the X direction
%
% This indexing scheme can be visualized as below
%      
%                     (i-C)   (i+R*C)
%                       |       /
%                       |     /
%                       |   /
%                       | / 
%       (i-1) <-------- i --------> (i+1)
%                     / |
%                   /   |
%                 /     |
%               /       |
%           (i-R*C)   (i+C)
%
% The volume looked at can be visualized by a cross-section
%
%        ____W____ 
%       |         |         y      z 
%       |         |         |    /
%       |         | H       |  /
%       |_________|         |/_____ x
%
%


function DF_SS_VK_3D
% System parameters
K1 = 5;                      % m/d
LENGTH = 1;                  % m
HEIGHT = 1;                  % m
WIDTH = 1;                   % m  
AREA = HEIGHT*WIDTH;         % m^2
Hin = 28;                    % m
Hout = 27.8;                 % m

% Model Parameters
R = 20;            % Nodes per Row
C = 20;            % Nodes per Column
L = 20;            % Nodes per Layer
N = R*C*L;          % Total number of nodes
dx = WIDTH/C;       % m
dy = HEIGHT/R;      % m 
dz = LENGTH/L;      % m

% Compute these in advance for simplicity
dx2 = dx*dx;
dy2 = dy*dy;
dz2 = dz*dz;

% Construct permeability field
K = zeros(N,1);
for i=1:C
    for j=1:R
        for k=1:L
            IC = R*C*(k-1)+C*(j-1)+i;
            K(IC) = K1;
        end
    end
end

%%%%%%%%%%% Start the calculations to determine the head field %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Allocate for M and b
M = zeros(N,N);     % Model matrix of discretized equations
b = zeros(N,1);     % RHS of model

% Populate b
for i=1:C
    for j=1:R
        LI = R*C*(1-1)+C*(j-1)+i;
        RI = R*C*(L-1)+C*(j-1)+i;
        b(LI) = -(K(LI)/(dz2) )*Hin;            % Inlet
        b(RI) = -(K(RI)/(dz2) )*Hout;           % Outlet
    end
end

% Populate M
for i=1:C
    for j=1:R
        for k=1:L
            % Prepare indices
            IC = R*C*(k-1)+C*(j-1)+i;
            PZ = IC + R*C; PY = IC + C; PX = IC + 1;
            MZ = IC - R*C; MY = IC - C; MX = IC - 1;
            
            % Calculate Effective Viscosities
            if (i>1), K_MX = K_Eff(K(IC),K(MX)); end     % Excludes Left 
            if (i<C), K_PX = K_Eff(K(IC),K(PX)); end     % Excludes Right 
            if (j>1), K_MY = K_Eff(K(IC),K(MY)); end     % Excludes Top 
            if (j<R), K_PY = K_Eff(K(IC),K(PY)); end     % Excludes Bottom
            if (k>1), K_MZ = K_Eff(K(IC),K(MZ)); end     % Excludes Front
            if (k<L), K_PZ = K_Eff(K(IC),K(PZ)); end     % Excludes Back
            
            % Populate M
            if (i>1), M(IC,MX) = K_MX/dx2; M(IC,IC) = M(IC,IC)-K_MX/dx2; end % Excludes Left 
            if (i<C), M(IC,PX) = K_PX/dx2; M(IC,IC) = M(IC,IC)-K_PX/dx2; end % Excludes Right 
            if (j>1), M(IC,MY) = K_MY/dy2; M(IC,IC) = M(IC,IC)-K_MY/dy2; end % Excludes Top
            if (j<R), M(IC,PY) = K_PY/dy2; M(IC,IC) = M(IC,IC)-K_PY/dy2; end % Excludes Bottom
            if (k>1), M(IC,MZ) = K_MZ/dz2; M(IC,IC) = M(IC,IC)-K_MZ/dz2; end % Excludes Front
            if (k<L), M(IC,PZ) = K_PZ/dz2; M(IC,IC) = M(IC,IC)-K_PZ/dz2; end % Excludes Back
            
            % Add Dirichlet front and back
            if (k==1), M(IC,IC) = M(IC,IC)-K(IC)/dz2; end % Front
            if (k==L), M(IC,IC) = M(IC,IC)-K(IC)/dz2; end % Back
        end
    end
end

H = linsolve(M,b);

visualPerp(H,R,C,dx,dy,2);
visualPar(H,R,C,L,dy,dz,2);

%%%%%%%%%%%% End the calculation of the head field %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

function K_Eff = K_Eff(K_1,K_2) 
    K_Eff = 2/(1/K_1 + 1/K_2); 
end

% Visualize Head Values in a plane perpendicular to flow
function visualPerp(H,R,C,dx,dy,k)
% Create axes
X = zeros(C,1);
    for i=1:C , X(i) = i*dx; end
Y = zeros(R,1);
    for j=1:R , Y(j) = j*dy; end

% Create viewable H matrix
H_View = zeros(R,C);
for i=1:C
    for j=1:R
        IC = R*C*(k-1)+C*(j-1)+i;
        H_View(j,i) = H(IC);
    end
end

figure('Name','View Perpendicular to Flow')
pcolor(X,Y,H_View);
ylabel(colorbar,'Pressure Head (m)');

end

% Visualize Head Values in a plane parallel to flow
function visualPar(H,R,C,L,dy,dz,i)
% Create axes
Z = zeros(L,1);
    for k=1:L , Z(k) = k*dz; end
Y = zeros(R,1);
    for j=1:R , Y(j) = j*dy; end

% Create viewable H matrix
H_View = zeros(R,L);
for k=1:L
    for j=1:R
        IC = R*C*(k-1)+C*(j-1)+i;
        H_View(j,k) = H(IC);
    end
end

figure('Name','View Parallel to Flow')
pcolor(Z,Y,H_View);
ylabel(colorbar,'Pressure Head (m)');

end