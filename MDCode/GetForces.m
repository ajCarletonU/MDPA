function GetForces(PhiCutoff,Epsilon,sigma)
global x y nAtoms Fx Fy Phi

% n = 1;

%1 LOOP for nAtoms
for i = 1:nAtoms
    Fx(i) = 0;
    Fy(i) = 0;
    Phi(i) = 0;
%     if i == 10
%         plot(x(i),y(i),'o','markers',48);
%         hold on
%     end


    %2 Loop (root atom )- for every single atom have to visit every other atom 
    % calc the forces and update the forces
    for j = 1:nAtoms
        if i == j, continue; end
        dx = x(i) - x(j);
        dy = y(i) - y(j);
        r = sqrt(dx^2 + dy^2);

        %only calc the force if the atom is relativily close to the root
        %atom (r) 
        if r > PhiCutoff, continue, end

        [aPhi dPhidr] = LJPot(r, Epsilon, sigma);

        ang = atan2(dy, dx);
%         dx =  r*cos(ang)
%         dy =  r*sin(ang)

        %updating forces
        dFx = - dPhidr * cos(ang);
        dFy = - dPhidr * sin(ang);

%         if i == 10
%             plot(x(j),y(j),'ro','markers',48);
%         end

        %adding all the forces on a particular atom
        Phi(i) = Phi(i) + aPhi;
        Fx(i) = Fx(i) + dFx;
        Fy(i) = Fy(i) + dFy;
    end
end

hold off

end

