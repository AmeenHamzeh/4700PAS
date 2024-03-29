set(0, 'DefaultFigureWindowStyle', 'docked')
set(0, 'defaultaxesfontsize', 18)
set(0, 'defaultaxesfontname', 'Times New Roman')
set(0, 'DefaultLineLineWidth', 2);

nx = 50;
ny = 50;
alpha = 1;
V = zeros(nx,ny);
G = sparse(nx*ny, nx*ny);

Inclusion = 0;

% for i = 1:nx
%     if i == 1 || i == nx
%         G(i,:) = 0;
%         G(i, i) = 1;
%     end 
% end
% 
% % Left and right boundary
% for j = 1:ny
%     if j == 1 || j == ny
%         G(j,:) = 0;
%         G(j, j) = 1;
%     end
% end

% for i = 1:nx
%     G(i, i) = 1; % Top row
%     G(i + (ny-1)*nx, i + (ny-1)*nx) = 1; % Bottom row
% end
% 
% % Left and right boundary
% for j = 1:ny
%     G((j-1)*nx+1, (j-1)*nx+1) = 1; % Left column
%     G(j*nx, j*nx) = 1; % Right column
% end


% Bulk nodes
% for i = 2:nx-1
%     for j = 2:ny-1
%         m = (j-1)*nx + i; % Node index in the single vector representation
% 
%         % Applying the Finite Difference equation
%         % Assuming the FD equation is Laplace's equation: d^2u/dx^2 + d^2u/dy^2 = 0
%         % with a 5-point stencil approximation
% 
%         G(m, m) = -4; % The diagonal element
%         G(m, m-1) = 1; % The left neighbor
%         G(m, m+1) = 1; % The right neighbor
%         G(m, m-nx) = 1; % The upper neighbor
%         G(m, m+nx) = 1; % The lower neighbor
%     end
% end
for i = 1:nx
    for j = 1:ny
        m = (i-1)*ny + j;
        if i == 1 || i == nx || j == 1 || j == ny  % Boundary nodes
            G(m, m) = 1;
            %G(m, :) = 0;
        else  % Internal nodes
              if i > 10 && i < 20 && j > 10 && j < 20  % Modified region
                  G(m, m) = -2*alpha;
              else
                  G(m,m) = -4;
              end
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;
            
            G(m,nxm) = 1;
            G(m,nxp) = 1;
            G(m,nym) = 1;
            G(m,nyp) = 1;
        
    end
    end
end

figure('name', 'Matrix')
spy(G)

nmodes = 20;
[E,D] = eigs(G,nmodes,'SM');

figure('name', 'EigenValues')
plot(diag(D),'*');

np = ceil(sqrt(nmodes));

figure('name', 'Modes')
for k = 1:nmodes
    M = E(:,k);
    for i = 1:nx
        for j = 1:ny
            n = j + (i-1)*nx;
            V(i,j) = M(n);
        end
    
    subplot(np,np,k), surf(V,'Linestyle','none')
    title(['EV =' num2str(D(k,k))])
    end
end