clear all;
format long;
% Flags
flagDBC = 1;        % flag to apply Dirichlet BC
flagPlotDoping = 1; % flag to plot doping profile in the beginning
flagPlotAll = 1;    % flag to plot the reults after reaching the convergence 
flagAnSol = 1;      % flag to calculate analytical solutions
nIter = 0;        	% index of Newton iteration

% Material Constants and Simulation Environment
eps0 = 8.85E-12;    % vacuum dielectric constant (F/m)
epsSi = 11.7;       % the relative dielectric constant of Si
eps = epsSi * eps0; % the dielectric constant (F/m) of Si
q0 = 1.602E-19;     % elementray charge (C)
ni = 1.5E10 * 1E6;  % the intrinsic carrier concentration (1/m3) of Si 
k0 = 8.617E-5;      % Boltzmann constant (eV/K)
T = 300;            % Temperature
VT = T*k0;          % thermal energy (eV)
Na = 1E16 * 1E6;    % p-doping concentration
Nd = 1E16 * 1E6;    % n-doping concentration

quasi_n = 0;        % quasi-Fermi level for electrons (prepared for future non-eqilibrium simulations)
quasi_p = 0;        % quasi-Fermi level for holes (prepared for future non-eqilibrium simulations)
rn = q0 * ni / VT * exp(-quasi_n/VT); % prefactor of electron concentration 
rp = q0 * ni / VT * exp( quasi_p/VT); % prefactor of hole concentration 

% Define the simulation domain 
xLen = 2E-6;        % the length of the PN diode (2 um)
yLen = 2E-6;        % the thicness of the PN diode (2 um)
NxElm = 10;         % number of elements along the x direction 
NyElm = 30;         % number of elements along the y direction 
NxNode = NxElm + 1; % number of nodes along the x direction 
NyNode = NyElm + 1; % number of nodes along the y direction 
dx = xLen/NxElm;    % length of each element along x (uniform mesh)
dy = yLen/NyElm;    % length of each element along y (uniform mesh)
totNElm = 0;        % initialize the total number of elements
totNode = NxNode*NyNode; % total number of nodes


% Initialization of [K], [b], and [V]
K = zeros(totNode,totNode); % matrix on the left-hand  side of Eq. 7
B = zeros(totNode,1);       % vector on the right-hand side of Eq. 7
V = zeros(totNode,1);       % Dirichlet BC
D = zeros(totNElm);         % source term related to the ionized impurities
Phi_old = zeros(totNode,1); % old potentials
Phi_new = zeros(totNode,1); % updated potentials
n = [ ];                    % store the nodal indices for each element
nd = [ ];                   % store the nodal indices on the Dirichlet BC
ND = 0;                     % initialized the number of nodes on the Dirichlet BC
idx2Dto1D = zeros(NxNode,NyNode); % initialize the table to map 2D indices onto 1D indices


% Initialize the nodal coordinates 
coord = struct('x',zeros(NxNode),'y', zeros(NyNode));  % use struct to store (x,y) coordinate information
r = 0.8; 		% setupthe non-uniform meshes along the y direction in the p-region (eq. 14-15)
h = (r - 1)/(r^(NyNode/2-1) - 1); 
ytmp1 = [0, h*cumsum(r.^(0:(NyNode/2-2)))];

r = 1.2;		% setup the non-uniform meshes along the y direction in the n-region (eq. 14-15)
h = (r - 1)/(r^((NyNode/2+1-1) - 1)); 
ytmp2 = [1, 1+h*cumsum(r.^(0:(NyNode/2-2)))];
yall = [ytmp1 ytmp2 2] * 0.5 * yLen;

for i = 1:NxNode
    for j = 1:NyNode
        coord(i,j) = struct('x', dx * (i-1),'y', yall(j)); % assign (x,y) coordinates to the nodes
    end 
end


% Assemble the nodes into the elements
for i = 1:NxNode-1
    for j = 1:NyNode-1
    if mod(i+j,2) == 0  
	% Define two triangular elements in a rectangle where the slope of the bounary is positive 
        p1 = struct('i', i,'j', j);
        p2 = struct('i', i+1,'j', j);
        p3 = struct('i', i+1,'j', j+1);
        n = [ n; p1 p2 p3];
        p2 = p3;
        p3 = struct('i', i,'j', j+1);
        n = [ n; p1 p2 p3];
    else
	% Define two triangular elements in a rectangle where the slope of the bounary is negative 
        p1 = struct('i', i,'j', j);
        p2 = struct('i', i+1,'j', j);
        p3 = struct('i', i,'j', j+1);
        n = [ n; p1 p2 p3];
        p1 = p2;
        p2 = struct('i', i+1,'j', j+1);
        n = [ n; p1 p2 p3];
    end
    totNElm = totNElm + 2; 
    end
end


% Convert 2D (i,j) coordinates to 1D indicies   
for i = 1:NxNode
    for j = 1:NyNode
          idx2Dto1D(i,j) = i + (j-1) * NxNode; 
    end
end


% Initialize the Dirichlet BC
if flagDBC == 1 
    % Specify the nodes with Dirichlet BC (the top and bottom of the simulation domain)
    for i = 1:NxNode
        p1 = struct('i', i,'j', 1);
        p2 = struct('i', i,'j', NyNode);
        nd = [nd; p1; p2];
        ND = ND+2;
    end
	
    % Setup the potentials at the Dirichlet BC  
    for i = 1:ND
        p1 = idx2Dto1D(nd(i).i, nd(i).j); 
        if  nd(i).j == 1 % the Dirichlet BC of the p-region
            V(p1) = log((-Na/2 + sqrt(Na^2/4 + ni^2))/ni);
        else             % the Dirichlet BC of the n-region
            V(p1) = log((Nd/2 + sqrt(Nd^2/4 + ni^2))/ni); 
        end
        Phi_old(p1) = V(p1); % initialize the potentials based on the nearby Dirichlet BC

    end
end

% Setup the soruce term related to the ionized impurities 
for e = 1:totNElm
        D(e) = q0 * Nd / VT;  
end
for e = 1:totNElm
    if n(e,1).j <= ceil(NyNode/2) && n(e,2).j <= ceil(NyNode/2) && n(e,3).j <= ceil(NyNode/2) 
        D(e) = -q0 * Na / VT;  
    end
end

%  Setup the initial guess of potentials 
%  The values are equal to the nearby Dirichlet BC (fig. 7)
for i = 1:totNode
    j = ceil(i/NxNode);    
    if  j <= ceil(NyNode/2)
        Phi_old(i) = log((-Na/2 + sqrt(Na^2/4 + ni^2))/ni);  % for p-region
    else
        Phi_old(i) = log((Nd/2 + sqrt(Nd^2/4 + ni^2))/ni);   % for n-region
    end
end



while(1)
K = zeros(totNode,totNode);  
J = zeros(totNode,totNode);
B = zeros(totNode,1);

% Assemble matrices for eqs. 7  ([K][psi] = [B])
for e = 1: totNElm
 
	a = n(e,1); 
    b = n(e,2);
    c = n(e,3);
    ae = zeros(3);
    be = zeros(3);
    ce = zeros(3);
    
    ae(1) = coord(b.i,b.j).x*coord(c.i,c.j).y - coord(c.i,c.j).x*coord(b.i,b.j).y;
    ae(2) = coord(c.i,c.j).x*coord(a.i,a.j).y - coord(a.i,a.j).x*coord(c.i,c.j).y; 
    ae(3) = coord(a.i,a.j).x*coord(b.i,b.j).y - coord(b.i,b.j).x*coord(a.i,a.j).y;
    be(1) = coord(b.i,b.j).y - coord(c.i,c.j).y;
    be(2) = coord(c.i,c.j).y - coord(a.i,a.j).y;
    be(3) = coord(a.i,a.j).y - coord(b.i,b.j).y;
    ce(1) = coord(c.i,c.j).x - coord(b.i,b.j).x;
    ce(2) = coord(a.i,a.j).x - coord(c.i,c.j).x;
    ce(3) = coord(b.i,b.j).x - coord(a.i,a.j).x;
    deltae = 0.5 * (be(1)*ce(2)-be(2)*ce(1)); % the area of the element
    
	
	% Assemble the K matrix element-by-element (LHS of eq. 7)
    Ke = zeros(3,3); % sub-K matrix for each element
    for i = 1:3
        for j = 1:3
            Ke(i,j) = eps * (be(i)*be(j) + ce(i)*ce(j))/(4*deltae); 
            p1 = idx2Dto1D(n(e,i).i, n(e,i).j);  % find the corresponding row in the K matrix
            p2 = idx2Dto1D(n(e,j).i, n(e,j).j);  % find the corresponding column in the K matrix
            K(p1,p2) = K(p1,p2) + Ke(i,j);
        end 
    end
	
	% Assemble the B vector element-by-element 
    for i = 1:3 % the 3 nodes (l,m,n) in the element 'e' 
            i1 = i; 
            i2 = mod(i1,3) + 1; 
            i3 = mod(i2,3) + 1; 
            
            p1 = idx2Dto1D( n(e,i1).i, n(e,i1).j );  % node l of eq. 7 
            p2 = idx2Dto1D( n(e,i2).i, n(e,i2).j );  % node m of eq. 7 
            p3 = idx2Dto1D( n(e,i3).i, n(e,i3).j );  % node n of eq. 7 
            
			% numerical integral for the electron part (the first term on the RHS) of eq. 7
            fun = @(x,y) x .* exp(-x*Phi_old(p1)-y*Phi_old(p2)-(1-x-y)*Phi_old(p3));
            ymax = @(x) 1 - x;
            B(p1) = B(p1) + rp*2*deltae*integral2(fun,0,1,0,ymax);
			
			% numerical integral for the hole part (the second term on the RHS) of eq. 7
            fun = @(x,y) x .* exp(x*Phi_old(p1)+y*Phi_old(p2)+(1-x-y)*Phi_old(p3));
            ymax = @(x) 1 - x;                  
            B(p1) = B(p1) - rn * 2*deltae*integral2(fun,0,1,0,ymax);
			
			% the contribution from the ionized impurities (the last term on the RHS of eq. 7)
			B(p1) = B(p1) + D(e)*deltae/3; 
    end
end

% Assign the Dirichlet BC by modifying K and B matrices
if flagDBC == 1 
    for d = 1:ND
        pD = idx2Dto1D( nd(d).i, nd(d).j); 
        for i = 1:NxNode
            for j = 1:NyNode
                p1 = idx2Dto1D(i,j); 
                if p1 == pD
                    K(p1,p1) = 1;
                    B(p1) = V(p1); 
                else
                    B(p1) = B(p1) - K(p1,pD) * V(pD); 
                    K(pD,p1) = 0;
                    K(p1,pD) = 0;
                end
            end
        end
    end
end


% Assemble the Jabobian matrix for eqs. 12 ([J])
for e = 1: totNElm
    a = n(e,1);
    b = n(e,2);
    c = n(e,3);
    ae = zeros(3);
    be = zeros(3);
    ce = zeros(3);
    ae(1) = coord(b.i,b.j).x*coord(c.i,c.j).y - coord(c.i,c.j).x*coord(b.i,b.j).y;
    ae(2) = coord(c.i,c.j).x*coord(a.i,a.j).y - coord(a.i,a.j).x*coord(c.i,c.j).y; 
    ae(3) = coord(a.i,a.j).x*coord(b.i,b.j).y - coord(b.i,b.j).x*coord(a.i,a.j).y;
    be(1) = coord(b.i,b.j).y - coord(c.i,c.j).y;
    be(2) = coord(c.i,c.j).y - coord(a.i,a.j).y;
    be(3) = coord(a.i,a.j).y - coord(b.i,b.j).y;
    ce(1) = coord(c.i,c.j).x - coord(b.i,b.j).x;
    ce(2) = coord(a.i,a.j).x - coord(c.i,c.j).x;
    ce(3) = coord(b.i,b.j).x - coord(a.i,a.j).x;
    deltae = 0.5 * (be(1)*ce(2)-be(2)*ce(1));
    
	% Assemble the first term on the RHS of eq. 12
    Je = zeros(3,3);
    for i = 1:3
        for j = 1:3
            Je(i,j) = eps * (be(i)*be(j) + ce(i)*ce(j))/(4*deltae); 
            p1 = idx2Dto1D(n(e,i).i, n(e,i).j);  % from 2D nodes to 1D node number
            p2 = idx2Dto1D(n(e,j).i, n(e,j).j);
            J(p1,p2) = J(p1,p2) + Je(i,j);
        end 
    end
    
    % Assemble the first term on the RHS of eq. 12
    Je = zeros(3,3);
    for i = 1:3
        for j = 1:3
            if i == j % Diagonal 
                 i1 = i; 
                 i2 = mod(i1,3) + 1; 
                 i3 = mod(i2,3) + 1; 
            
                 p1 = idx2Dto1D( n(e,i1).i, n(e,i1).j );  
                 p2 = idx2Dto1D( n(e,i2).i, n(e,i2).j );  
                 p3 = idx2Dto1D( n(e,i3).i, n(e,i3).j );  
                 
                 fun = @(x,y) x.* x .* exp(-x*Phi_old(p1)-y*Phi_old(p2)-(1-x-y)*Phi_old(p3));
                 ymax = @(x) 1 - x;
                 Je(i,j) = rp * 2*deltae*integral2(fun,0,1,0,ymax);
                 
                 fun = @(x,y) x.* x .* exp(x*Phi_old(p1)+y*Phi_old(p2)+(1-x-y)*Phi_old(p3));
                 ymax = @(x) 1 - x;
                 Je(i,j) = Je(i,j) + rn * 2*deltae*integral2(fun,0,1,0,ymax);
                 
            else    % Off-diagonal 
                 i1 = i; 
                 i2 = j; 
                 i3 = 6-i-j; 
            
                 p1 = idx2Dto1D( n(e,i1).i, n(e,i1).j );  
                 p2 = idx2Dto1D( n(e,i2).i, n(e,i2).j );  
                 p3 = idx2Dto1D( n(e,i3).i, n(e,i3).j );  
                 
                 fun = @(x,y) x.* y .* exp(-x*Phi_old(p1)-y*Phi_old(p2)-(1-x-y)*Phi_old(p3));
                 ymax = @(x) 1 - x;
                 Je(i,j) = rp * 2*deltae*integral2(fun,0,1,0,ymax);        
                 
                 fun = @(x,y) x.* y .* exp(x*Phi_old(p1)+y*Phi_old(p2)+(1-x-y)*Phi_old(p3));
                 ymax = @(x) 1 - x;
                 Je(i,j) = Je(i,j) + rn * 2*deltae*integral2(fun,0,1,0,ymax);     
            end

            p1 = idx2Dto1D(n(e,i).i, n(e,i).j);  % find the corresponding row in the J matrix
            p2 = idx2Dto1D(n(e,j).i, n(e,j).j);  % find the corresponding column in the J matrix
            J(p1,p2) = J(p1,p2) + Je(i,j);
        end 
    end
end


% Newton-Raphson Iterations
F = K * Phi_old - B;                   % F associates to the delta x of eq. 9
CorrectionMatrix = J^-1 * F;      
if flagDBC == 1 			           % Nodes on the Dirichlet BCs are not updated
    for d = 1:ND
        pD = idx2Dto1D( nd(d).i, nd(d).j);
        CorrectionMatrix(pD) = 0;  
    end
end
Phi_new = Phi_old - CorrectionMatrix;  % Update the new potentials associated to eq. 10
maxErr = max(abs(CorrectionMatrix));   % Find Out the maximum error
disp('MaxError and nIter');
disp(maxErr*VT);
disp(nIter);
if maxErr < 1E-6                       % Convergence:condition: the maximum error < 1 miron in unit of VT
    break;
end
Phi_old = Phi_new;                     % Assign new potentials to old potentials and proceed to the next iteration
nIter = nIter + 1;
end

% Plot the initial doping profiles
if flagPlotDoping == 1
    figure(1);
    title('Initial Doping Profile');
    hold on;
    for e = 1: totNElm
        if D(e) > 0
            plot(polyshape([coord(n(e,1).i,n(e,1).j).x coord(n(e,2).i,n(e,2).j).x ...
            coord(n(e,3).i,n(e,3).j).x], [coord(n(e,1).i,n(e,1).j).y coord(n(e,2).i,n(e,2).j).y ...
            coord(n(e,3).i,n(e,3).j).y]), 'FaceAlpha', 1, 'FaceColor', 'red' );
        else
            plot(polyshape([coord(n(e,1).i,n(e,1).j).x coord(n(e,2).i,n(e,2).j).x ...
            coord(n(e,3).i,n(e,3).j).x], [coord(n(e,1).i,n(e,1).j).y coord(n(e,2).i,n(e,2).j).y ...
            coord(n(e,3).i,n(e,3).j).y]), 'FaceAlpha', 1, 'FaceColor', 'blue' );
        end
    end
    hold off;
end 

%Plot important results
if flagPlotAll == 1
    
    Phi_old2D = zeros(NxNode,NyNode);
    Phi_new2D = zeros(NxNode,NyNode);
    
    % Map 1D array into 2D array
    for index1D = 1:totNode
        j = ceil(index1D/NxNode); 
        i = index1D - (j-1)*NxNode;
        Phi_new2D(i,j) = Phi_new(index1D); 
        Phi_old2D(i,j) = Phi_old(index1D); 
    end

	% Plot the triangular meshes
    figure(2);
    title('Non-Uniform Mesh Generation');
    hold on;
    for e = 1: totNElm
        plot(polyshape([coord(n(e,1).i,n(e,1).j).x coord(n(e,2).i,n(e,2).j).x ...
            coord(n(e,3).i,n(e,3).j).x], [coord(n(e,1).i,n(e,1).j).y coord(n(e,2).i,n(e,2).j).y ...
            coord(n(e,3).i,n(e,3).j).y]), 'FaceAlpha', 0);
    end
    axis([0 2E-6 0 2E-6]);
    xlabel('x (m)');
    ylabel('y (m)');
    hold off;
    
	% Plot the converged potentials 
    figure(3);
    title('Phi_{new}');
    hold on;
    contourf( [0:dx:xLen], yall ,Phi_new2D'*VT);
    for e = 1: totNElm
         plot(polyshape([coord(n(e,1).i,n(e,1).j).x coord(n(e,2).i,n(e,2).j).x ...
            coord(n(e,3).i,n(e,3).j).x], [coord(n(e,1).i,n(e,1).j).y coord(n(e,2).i,n(e,2).j).y ...
            coord(n(e,3).i,n(e,3).j).y]), 'FaceAlpha', 0);
    end
    hold off;
	
end 

%Calculate the anaytical solutions for the comparison
if flagAnSol == 1
    PhiBuiltin = VT * log(Nd*Na/ni/ni);
    Wn = sqrt(2*eps*PhiBuiltin/q0 * Na/Nd * (Na+Nd)^-1);
    Wp = sqrt(2*eps*PhiBuiltin/q0 * Nd/Na * (Na+Nd)^-1);
    W = sqrt(2*eps*PhiBuiltin/q0 * (1/Na+1/Nd));
    tmpxp = -Wp:1E-9:0;
    tmpxn = 1E-9:1E-9:Wn;
    Phip = q0*Na/eps * 1/2 * (Wp + tmpxp).^2;
    Phin = max(Phip) -  q0*Nd/eps * 1/2 * tmpxn .* (tmpxn - 2*Wn);
    xNP = [tmpxp tmpxn];
    PhinNP = [Phip Phin];
    
end
