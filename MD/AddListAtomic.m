function AddListAtomic(X0, Y0, VX0, VY0, Type, InitDist, Temp)
global C
global x y AtomSpacing
global nAtoms
global AtomType Vx Vy Mass0 Mass1 Mass2

xp = X0;
yp = Y0;
numAtoms = length(X0);

x(nAtoms + 1:nAtoms+numAtoms) = xp + (rand(1, numAtoms) - 0.5) * AtomSpacing * InitDist;
y(nAtoms + 1:nAtoms+numAtoms) = yp + (rand(1, numAtoms) - 0.5) * AtomSpacing * InitDist;

% iTypes = Types == 1;
% Mass = ones(1,numAtoms)*Mass0;
% Mass(iTypes) = Mass1;

if Type == 2
    Mass = Mass2;
elseif Type == 0
    Mass = Mass0;
else
    Mass = Mass1;
end

AtomType(nAtoms + 1:nAtoms + numAtoms) = Type;


if Temp == 0
    Vx(nAtoms + 1:nAtoms + numAtoms) = 0;
    Vy(nAtoms + 1:nAtoms + numAtoms) = 0;
else
    std0 = sqrt(C.kb * Temp ./ Mass);

    Vx(nAtoms + 1:nAtoms + numAtoms) = std0 .* randn(1, numAtoms);
    Vy(nAtoms + 1:nAtoms + numAtoms) = std0 .* randn(1, numAtoms);
end

Vx(nAtoms + 1:nAtoms + numAtoms) = Vx(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vx(nAtoms + 1:nAtoms + numAtoms)) + VX0;
Vy(nAtoms + 1:nAtoms + numAtoms) = Vy(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vy(nAtoms + 1:nAtoms + numAtoms)) + VY0;

nAtoms = nAtoms + numAtoms;

end
