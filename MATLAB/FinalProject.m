clc; clear;
format long

bodies = buildSystem();
bodies = [bodies(1),bodies(2),bodies(3),bodies(4)];
n = 1e3;
dt = 60*60*24*365/n;
n = n*1;
lim = 0;
method = 'FWEuler';
graph = 1;
reverse = 0;

nrg = subplot(2,2,[3,4]);
E0 = getEnergy(bodies);
plot(nrg,0,0);
if strcmp(method,'FWEuler')
    title("Energy vs Time (Forward Euler)");
end
if strcmp(method,'BWEuler')
    title("Energy vs Time (Backward Euler)");
end
if strcmp(method,'EulerCromer')
    title("Energy vs Time (Euler-Cromer)");
end
if strcmp(method,'LeapFrog')
    title("Energy vs Time (Leap-Frog)");
end
if strcmp(method,'RKN89')
    title("\textbf{Energy vs Time (RKN8(9)-13-$\dot x$)}",'Interpreter', 'latex');
end
xlabel("time (s)");
ylabel("REL \DeltaE");
xlim(nrg,[0,n*dt]);
axis(nrg,'auto y');
hold on

sys = subplot(2,2,[1,2]);
plot3(0,0,0,'w');
if strcmp(method,'FWEuler')
    title("Forward Euler");
end
if strcmp(method,'BWEuler')
    title("Backward Euler");
end
if strcmp(method,'EulerCromer')
    title("Euler-Cromer");
end
if strcmp(method,'LeapFrog')
    title("Leap-Frog");
end
if strcmp(method,'RKN89')
    title("\textbf{RKN8(9)-13-$\dot x$}",'Interpreter', 'latex');
end
axis(sys,[-5e11,5e11,-5e11,5e11,-5e11,5e11]);
%axis([-1e13,1e13,-1e13,1e13,-1e13,1e13]);
hold on
for i=1:length(bodies)
    r = bodies(i).getPosition();
    %plot3(sys,r(1),r(2),r(3),'Color',bodies(i).getColor(),'Marker','.');
    plot3(sys,r(1),r(2),r(3),'r*');
end
for k = 1:n
    if strcmp(method,'FWEuler')
        stepFWEuler(bodies, dt);
    end
    if strcmp(method,'BWEuler')
        stepBWEuler(bodies, dt);
    end
    if strcmp(method,'EulerCromer')
        stepEulerCromer(bodies,dt);
    end
    if strcmp(method,'LeapFrog')
        stepLeapFrog(bodies, dt);
    end
    if strcmp(method,'RKN89')
        stepRKN89(bodies, dt);
    end
    if (k >= floor(lim*n))
        if (graph)
            for i = 1:length(bodies)
                r = bodies(i).getPosition();
                plot3(sys,r(1),r(2),r(3),'Color',bodies(i).getColor(),'Marker','.');
            end
        end
        lim = lim + 0.0005;
        
        E = getEnergy(bodies);
        plot(nrg,k*dt,abs((E0 - E))/E,'r.');
    end
    
        %Reverse Time halfway
     if k == n*0.5-1 && reverse
         for i=1:length(bodies)
             setVelocity(bodies(i),-bodies(i).getVelocity());
         end
     end
end
abs((E0 - E))/E
