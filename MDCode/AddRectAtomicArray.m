function AddRectAtomicArray(LAtoms, WAtoms, X0, Y0, VX0, VY0, InitDist, Temp, Type)
global C
global x y AtomSpacing
global nAtoms
global AtomType Vx Vy Mass0 Mass1

if Type == 0
    Mass = Mass0;
else
    Mass = Mass1;
end

%Sets initial position of the ATOMS
%------------------------------------
L = (LAtoms - 1) * AtomSpacing;
W = (WAtoms - 1) * AtomSpacing;

numAtoms = LAtoms * WAtoms;

xp(1, :) = linspace(0, L, LAtoms);          %creates the position of the atoms
yp(1, :) = linspace(0, W, WAtoms);

x(nAtoms + 1:nAtoms+LAtoms) = xp-L/2;
y(nAtoms + 1:nAtoms+LAtoms) = yp(1)-W/2;

for i = 1:WAtoms-1
    x(nAtoms + i * LAtoms + 1:nAtoms + (i + 1) * LAtoms) = xp - L / 2;
    y(nAtoms + i * LAtoms + 1:nAtoms + (i + 1) * LAtoms) = yp(i + 1) - W / 2;
end


%Disturbing the positions of the ATOMS
%- using rand() a number generator
%- by taking the original initial positions (on a prefect Rectang matrix) 
%  and adding a small random position to the atoms
%------------------------------------
x(nAtoms + 1:nAtoms + numAtoms) = x(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms) - 0.5) * AtomSpacing * InitDist + X0;
y(nAtoms + 1:nAtoms + numAtoms) = y(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms) - 0.5) * AtomSpacing * InitDist + Y0;

AtomType(nAtoms + 1:nAtoms + numAtoms) = Type;


%Calc thermal velocities
%------------------------------------
if Temp == 0
    Vx(nAtoms + 1:nAtoms + numAtoms) = 0;  %velocity is zero if the Temp =o
    Vy(nAtoms + 1:nAtoms + numAtoms) = 0;
else
    std0 = sqrt(C.kb * Temp / Mass);       %calc the distribution

    %using randn() - is a normal distrubion function 
    Vx(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
    Vy(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
end


%Group velocities
%------------------------------------
%adding VX_knot and VY_knot to all the atoms
Vx(nAtoms + 1:nAtoms + numAtoms) = Vx(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vx(nAtoms + 1:nAtoms + numAtoms)) + VX0;
Vy(nAtoms + 1:nAtoms + numAtoms) = Vy(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vy(nAtoms + 1:nAtoms + numAtoms)) + VY0;

nAtoms = nAtoms + numAtoms;

end
