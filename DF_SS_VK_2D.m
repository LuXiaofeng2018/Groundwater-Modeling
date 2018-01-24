% This is a first attempt at calibrating a model for a system with
% multiple different permeabilities. The system will have one media on the
% left half of the system, a different media on the right half. The
% pressure drops across both sections should be determined analytically
% using Darcy's law and an effective permeability based on the harmonic
% mean.
%
% Rows go in the x-direction, columns in the y-direction.
% The model's indices are numbered from left to right, top to bottom such
% that, for node i in x and j in y, the center node is C*(j-1)+i. The
% indices which are used are labeled below
%
%       IC = C*(j-1)+i;     IC is the index for the center node
%       PY = IC + C;        PY is the index for +1 in the Y direction
%       PX = IC + 1;        PX is the index for +1 in the X direction
%       MY = IC - C;        MY is the index for -1 in the Y direction
%       MX = IC - 1;        MX is the index for -1 in the X direction
%
% This indexing scheme can be visualized as below
%                   
%                 (i-C)
%                   |
%                   |
%       (i-1) <---- i ----> (i+1)
%                   |
%                   |
%                 (i+C)
%
% v 0.1, first real go.

% Darcy flow, steady-state, variable K, 2D
function DF_SS_VK_2D
% System parameters
K1 = 5;                         % m/d
K2 = 5;                       % m/d
LENGTH = 1;                     % m
HEIGHT = 1;                     % m
WIDTH = 1;                      % m  
AREA = HEIGHT*WIDTH;            % m^2
Hin = 28;                       % m
Hout = 27.8;                    % m

% Model Parameters
R = 20;            % Nodes per Row
C = 20;            % Nodes per Column
N = R*C;            % Total number of nodes
dx = LENGTH/C;      % m
dy = HEIGHT/R;      % m  

% Compute these in advance for simplicity
dx2 = dx*dx;
dy2 = dy*dy;

% Construct permeability vector
K = Permeability(N,R,C,K1,K2);

% Compute the head values
H = Calculate(N,R,C,dx2,dy2,K,Hin,Hout);

% Visualize the data with homemade function
Visualize(H,R,C,dx,dy)

% Compare to analytical solution
Error_Check(K1,K2,K,AREA,Hin,Hout,H,LENGTH,dx,C)


end

% Generate field of permeability values
function K = Permeability(N,R,C,K1,K2)
K = zeros(N,1);
for i_K=1:C
    for j_K=1:R
        if (i_K<=C/2), K(C*(j_K-1)+i_K) = K1; end
        if (i_K>C/2),  K(C*(j_K-1)+i_K) = K2; end
    end
end


end

function K_Eff = K_Eff(K_1,K_2)
K_Eff = 2/(1/K_1 + 1/K_2);
end

% Begin the actual calculation 
function H = Calculate(N,R,C,dx2,dy2,K,Hin,Hout)

% Allocate for M and b
M = zeros(N,N);     % Model matrix of discretized equations
b = zeros(N,1);     % RHS of model

% Populate b
for j=1:R
    LI = 1+C*(j-1); RI = C+C*(j-1);     % Left and Right Index
    b(LI) = -(K(LI)/(dx2) )*Hin;            % Inlet
    b(RI) = -(K(RI)/(dx2) )*Hout;           % Outlet
end

% Apply boundary conditions
% Corner Nodes, starting with upper left
IC = 1; PY = IC + C; PX = IC + 1;

K_PX = K_Eff(K(IC),K(PX)); 
K_PY = K_Eff(K(IC),K(PY)); 
K_MX = K(IC);   

M(IC,IC) = -(K_PY/dy2 + K_PX/dx2 + K_MX/dx2 );
M(IC,PX) = K_PX/dx2;
M(IC,PY) = K_PY/dy2;

% Upper Right
IC = C; PY = IC + C; MX = IC - 1;

K_MX = K_Eff(K(IC),K(MX)); 
K_PY = K_Eff(K(IC),K(PY)); 
K_PX = K(IC);   

M(IC,IC) = -(K_PY/dy2 + K_PX/dx2 + K_MX/dx2 );
M(IC,MX) = K_MX/dx2;
M(IC,PY) = K_PY/dy2;

% Lower Left
IC = C*(R-1)+1; PX = IC + 1; MY = IC - C;

K_PX = K_Eff(K(IC),K(PX)); 
K_MY = K_Eff(K(IC),K(MY)); 
K_MX = K(IC);   

M(IC,IC) = -(K_MY/dy2 + K_PX/dx2 + K_MX/dx2 );
M(IC,PX) = K_PX/dx2;
M(IC,MY) = K_MY/dy2;

% Lower Right
IC = C*(R-1)+C; MY = IC - C; MX = IC - 1;

K_MX = K_Eff(K(IC),K(MX)); 
K_MY = K_Eff(K(IC),K(MY)); 
K_PX = K(IC);   

M(IC,IC) = -(K_MY/dy2 + K_PX/dx2 + K_MX/dx2 );
M(IC,MX) = K_MX/dx2;
M(IC,MY) = K_MY/dy2;

% Upper and Lower Middle Nodes
for i=2:(R-1)
    % Upper Nodes
    IC = i; PY = IC + C; PX = IC + 1; MX = IC - 1;
    
    K_MX = K_Eff(K(IC),K(MX)); 
    K_PY = K_Eff(K(IC),K(PY)); 
    K_PX = K_Eff(K(IC),K(PX));

    M(IC,IC) = -(K_PY/dy2 + K_PX/dx2 + K_MX/dx2 );
    M(IC,MX) = K_MX/dx2;
    M(IC,PY) = K_PY/dy2;
    M(IC,PX) = K_PX/dx2;
    
    % Lower Nodes
    IC = C*(R-1)+i; MY = IC - C; PX = IC + 1; MX = IC - 1;
    
    K_MX = K_Eff(K(IC),K(MX)); 
    K_MY = K_Eff(K(IC),K(MY)); 
    K_PX = K_Eff(K(IC),K(PX));
    
    M(IC,IC) = -(K_MY/dy2 + K_PX/dx2 + K_MX/dx2 );
    M(IC,MX) = K_MX/dx2;
    M(IC,MY) = K_MY/dy2;
    M(IC,PX) = K_PX/dx2;
end

% Left and Right Middle Nodes
for j=2:(R-1)
    % Left Nodes
    IC = C*(j-1)+1; PY = IC + C; MY = IC - C; PX = IC + 1; 
    
    K_MX = K(IC); 
    K_MY = K_Eff(K(IC),K(MY)); 
    K_PX = K_Eff(K(IC),K(PX));
    K_PY = K_Eff(K(IC),K(PY));
    
    M(IC,IC) = -(K_PY/dy2 + K_PX/dx2 + K_MX/dx2 + K_MY/dy2 );
    M(IC,MY) = K_MY/dy2;
    M(IC,PY) = K_PY/dy2;
    M(IC,PX) = K_PX/dx2;
    
    % Right Nodes
    IC = C*(j-1)+C; MY = IC - C; PY = IC + C; MX = IC - 1;
    
    K_PX = K(IC); 
    K_MY = K_Eff(K(IC),K(MY)); 
    K_MX = K_Eff(K(IC),K(MX));
    K_PY = K_Eff(K(IC),K(PY));
    
    M(IC,IC) = -(K_PY/dy2 + K_PX/dx2 + K_MX/dx2 + K_MY/dy2 );
    M(IC,MX) = K_MX/dx2;
    M(IC,MY) = K_MY/dy2;
    M(IC,PY) = K_PY/dy2;
end
% end boundary conditions

% Populate M
for i = 2:(C-1)
    for j = 2:(R-1)
        IC = C*(j-1)+i; PY = IC + C; PX = IC + 1; MY = IC - C; MX = IC - 1;
        
        K_MX = K_Eff(K(IC),K(MX)); 
        K_MY = K_Eff(K(IC),K(MY)); 
        K_PX = K_Eff(K(IC),K(PX));
        K_PY = K_Eff(K(IC),K(PY));
        
        M(IC,PY) = K_PY/dy2;
        M(IC,PX) = K_PX/dx2;
        M(IC,IC) = -(K_PY/dy2 + K_PX/dx2 + K_MX/dx2 + K_MY/dy2 );
        M(IC,MX) = K_MX/dx2;
        M(IC,MY) = K_MY/dy2;
    end
end

H = linsolve(M,b);

end

% Visualize Head Values
function Visualize(H,R,C,dx,dy)
% Create axes
X = zeros(C,1);
    for i=1:C , X(i) = i*dx; end
Y = zeros(R,1);
    for j=1:R , Y(j) = j*dy; end

% Create viewable H matrix
H_View = zeros(R,C);
for i=1:C
    for j=1:R
        H_View(j,i) = H(C*(j-1)+i);
    end
end

figure
pcolor(X,Y,H_View);
ylabel(colorbar,'Pressure Head (m)');

end

% Check for errors
function Error_Check(K1,K2,K,AREA,Hin,Hout,H,LENGTH,dx,C)
Kf = K_Eff(K1,K2);
A = AREA; L = LENGTH+dx;
Q = Kf*A*(Hin-Hout)/L;

H_Anal = Hin-Q*dx/(A*K1);
max_error = abs(H_Anal-H(1));
for i=2:C
    L = dx; Kf = K_Eff(K(i-1),K(i));
    H_Anal = H_Anal-Q*L/(A*Kf);
    error = abs(H_Anal-H(i));
    if (error>max_error), max_error=error; end
end

display(max_error);
end







