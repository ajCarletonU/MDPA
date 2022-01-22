%PA2 - ELEC4700 winter 2022
%Alina Jacobson (101055071) 

%Part2 choose to do para (a) Add a new geometry to the Add*Atomic*.m.
%elliptical
%baseline code used from AddCircAtomicArray.m
%added New function
%changes made to the postion of the atoms
%changes made to the equation of the cicle x^2 +y^2--> horizontal elliplical

function AddEllipticalAtomicArray(RADatomsA,RADatomsB,X0,Y0,VX0,VY0,InitDist,Temp,Type)

% initialized global variables
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
L = (2*RADatomsA - 1) * AtomSpacing;    % Lenght of the horizontal elliptical on the x-axis
W = (2*RADatomsB - 1) * AtomSpacing;    % width of the elliptical on the y-axis

xp(1, :) = linspace(-L/2, L/2, 2*RADatomsA);   %creates the position of the atoms
yp(1, :) = linspace(-W/2, W/2, 2*RADatomsB);


numAtoms = 0;           %set number of atoms at initial to be zero

%create elliptical
%loop through a range lenght
for i = 1:2*RADatomsA
    
    %loop through range height
    for j = 1:2*RADatomsB
      
        
        %updated equation    
        %horizontal elliplical equation = x^2/r^2 + y^2/r^2 = 1
        if(xp(i)^2/(RADatomsA*AtomSpacing)^2 + yp(j)^2/(RADatomsB*AtomSpacing)^2)<=1
            
            numAtoms = numAtoms+1;
            x(nAtoms + numAtoms) = xp(i);
            y(nAtoms + numAtoms) = yp(j);
        else
            i
            j
        end
    end
end



%Disturbing the positions of the ATOMS in the elliptical region
%- using rand() a number generator
%- by taking the original initial positions (on a prefect matrix) 
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

