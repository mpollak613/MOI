%* Copyright 2020 Michael Pollak.
%*
%* Use of this source code is governed by an MIT-style
%* licence that can be found in the LICENSE file.
%*

classdef Body < handle
    properties (Access = private)
       name string
       pos(1,3) double {mustBeReal, mustBeFinite}
       vel(1,3) double {mustBeReal, mustBeFinite}
       mass (1,1) double {mustBeReal, mustBeNonnegative, mustBeFinite}
       rad (1,1) double {mustBeReal, mustBeNonnegative, mustBeFinite}
       color (1,3) double
    end
    methods (Access = public)
        function obj = Body(name, pos, vel, mass, rad, col)
            if nargin == 6
                obj.name = name;
                obj.pos = pos;
                obj.vel = vel;
                obj.mass = mass;
                obj.rad = rad;
                obj.color = col;
            end
        end
        function [r] = getPosition(obj)
            r = obj.pos;
        end
        function [v] = getVelocity(obj)
           v = obj.vel;
        end
        function [c] = getColor(obj)
            c = obj.color;
        end
        function [n] = getName(obj)
            n = obj.name;
        end
        function [m] = getMass(obj)
            m = obj.mass;
        end
        function setVelocity(obj, nvel)
            obj.vel = nvel;
        end
        function newsys = copy(sys, delt)
            for i=1:length(sys)
                newsys(i) = Body('tmp',sys(i).pos + delt*sys(i).vel, ...
                    sys(i).vel,sys(i).mass,sys(i).rad,sys(i).color); %#ok<AGROW>
            end
        end
        function newsys = copy2(sys, npos)
            for i=1:length(sys)
                newsys(i) = Body('tmp',sys(i).pos + npos(i,:), ...
                    sys(i).vel,sys(i).mass,sys(i).rad,sys(i).color); %#ok<AGROW>
            end
        end
        function a = getForce(obj, bodies)
            G = 6.6742e-11; % m^3 * kg^-1 * s^-2
            a = 0;
            for i = setdiff(1:length(bodies),find(bodies==obj))
                a = a + G * bodies(i).mass * (bodies(i).pos - obj.pos)...
                    / norm(bodies(i).pos - obj.pos)^3;
            end
        end
        function E = getEnergy(bodies)
            G = 6.6742e-11; % m^3 * kg^-1 * s^-2
            K = 0;
            U = 0;
            for i = 1:length(bodies)
                for j = 1:i-1
                    K = K - G*bodies(j).mass*bodies(i).mass ...
                    /norm(bodies(i).pos - bodies(j).pos);
                end
                U = U + norm(bodies(i).mass*bodies(i).vel)^2 ...
                    /(2*bodies(i).mass);
            end
            E = K + U;
        end
        function stepFWEuler(bodies, delt)
            dr = zeros(length(bodies), 3);
            dv = zeros(length(bodies), 3);
            
            % FW Euler
            for i = 1:length(bodies)
                a = getForce(bodies(i), bodies);
                dr(i,:) = bodies(i).vel*delt;
                dv(i,:) = a*delt;
            end
            
            % update stats
            for i = 1:length(bodies)
                bodies(i).pos = bodies(i).pos + dr(i,:);
                bodies(i).vel = bodies(i).vel + dv(i,:);
            end
        end
        function stepBWEuler(bodies, delt)
            dr = zeros(length(bodies), 3);
            dv = zeros(length(bodies), 3);
            
            % BW Euler
            newsys = copy(bodies,delt);
            for i = 1:length(bodies)
                a = getForce(newsys(i), newsys);
                dv(i,:) = bodies(i).vel + a*delt;
                dr(i,:) = bodies(i).pos + dv(i,:)*delt;
            end
            
            % update stats
            for i = 1:length(bodies)
                bodies(i).pos = dr(i,:);
                bodies(i).vel = dv(i,:);
            end
        end
        function stepEulerCromer(bodies,delt)
            dr = zeros(length(bodies), 3);
            dv = zeros(length(bodies), 3);
            
            % Euler-Cromer
            for i=1:length(bodies)
                dv(i,:) = bodies(i).vel + getForce(bodies(i),bodies)*delt;
                dr(i,:) = bodies(i).pos + dv(i,:)*delt;
            end
            
            % update stats
            for i = 1:length(bodies)
                bodies(i).pos = dr(i,:);
                bodies(i).vel = dv(i,:);
            end
        end
        function stepLeapFrog(bodies, delt)
            dr = zeros(length(bodies), 3);
            dv = zeros(length(bodies), 3);
            
            % leap-frog
            newsys = copy(bodies, delt);
            for i = 1:length(bodies)
                a = getForce(bodies(i), bodies);
                a_1 = getForce(newsys(i), newsys);
                dr(i,:) = bodies(i).vel*delt + 0.5*a*delt^2;
                dv(i,:) = 0.5*(a + a_1)*delt;
            end
            
            % update stats
            for i = 1:length(bodies)
                bodies(i).pos = bodies(i).pos + dr(i,:);
                bodies(i).vel = bodies(i).vel + dv(i,:);
            end
        end
        function stepRKN89(bodies,delt)
            % formula obtained from
            % https://ntrs.nasa.gov/api/citations/19730015887/downloads/19730015887.pdf
            sqrt15 = sqrt(15);
            
            alpha = zeros(1,13);
            alpha(1) = 1/3;
            alpha(2) = 2/3;
            alpha(3) = 1/2;
            alpha(4) = 1/3;
            alpha(5) = 1;
            alpha(6) = 1/9;
            alpha(7) = 1/2;
            alpha(8) = 1/3;
            alpha(9) = 2/3;
            alpha(10) = 1/10 * (5 - sqrt15);
            alpha(11) = 1/10 * (5 + sqrt15);
            alpha(12) = 1;
            alpha(13) = 1;
            
            gamm = zeros(13,13);
            gamm(1,1) = 1/18;
            gamm(2,1) = 2/27; gamm(2,2) = 4/27;
            gamm(3,1) = 7/128; gamm(3,2) = 5/64; gamm(3,3) = -1/128;
            gamm(4,1) = 89/3240; gamm(4,2) = 31/540; ...
                gamm(4,3) = 11/1080; gamm(4,4) = -16/405;
            gamm(5,1) = 11/120; gamm(5,2) = 0; gamm(5,3) = 9/40; ...
                gamm(5,4) = -4/15; gamm(5,5) = 9/20;
            gamm(6,1) = 33259/7085880; gamm(6,2) = 0; ...
                gamm(6,3) = 343/157464; gamm(6,4) = -4708/885735; ...
                gamm(6,5) = 1879/393660; gamm(6,6) = -139/885735;
            gamm(7,1) = 29/1920; gamm(7,2) = 0; gamm(7,3) = 0; ...
                gamm(7,4) = 0; gamm(7,5) = 99/2560; ...
                gamm(7,6) = 1/30720; gamm(7,7) = 729/10240;
            gamm(8,1) = 13/1215; gamm(8,2) = 0; gamm(8,3) = 0; ...
                gamm(8,4) = 0; gamm(8,5) = 1/144; gamm(8,6) = 1/77760; ...
                gamm(8,7) = 87/2240; gamm(8,8) = -8/8505;
            gamm(9,1) = 22/1215; gamm(9,2) = 0; gamm(9,3) = 0; ...
                gamm(9,4) = 0; gamm(9,5) = 0; gamm(9,6) = 1/4860; ...
                gamm(9,7) = 3/28; gamm(9,8) = 256/8505; gamm(9,9) = 1/15;
            gamm(10,1) = 1/420000*(7561 - 1454*sqrt15); gamm(10,2) = 0; ...
                gamm(10,3) = 0; gamm(10,4) = 0; ...
                gamm(10,5) = -9/800000*(1373 + 45*sqrt15); ...
                gamm(10,6) = 1/1344000*(397 - 145*sqrt15); ...
                gamm(10,7) = 729/78400000*(6997 - 1791*sqrt15); ...
                gamm(10,8) = 1/183750*(999 - 473*sqrt15); ...
                gamm(10,9) = 27/5600000*(19407 - 3865*sqrt15); ...
                gamm(10,10) = 297/700000*(78 - 19*sqrt15);
            gamm(11,1) = 1/840000*(12647 + 2413*sqrt15); ...
                gamm(11,2) = 0; gamm(11,3) = 0; gamm(11,4) = 0; ...
                gamm(11,5) = -9/800000*(1373 - 45*sqrt15); ...
                gamm(11,6) = -1/336000*(29 - 61*sqrt15); ...
                gamm(11,7) = 729/19600000*(14743 + 3789*sqrt15); ...
                gamm(11,8) = 1/183750*(999 + 143*sqrt15); ...
                gamm(11,9) = 27/5600000*(20157 + 4315*sqrt15); ...
                gamm(11,10) = 27/1400000*(1641 + 463*sqrt15); ...
                gamm(11,11) = -1/56*(27 + 7*sqrt15);
            gamm(12,1) = 9/280 - 35/6561*5/2; gamm(12,2) = 0; ...
                gamm(12,3) = 0; gamm(12,4) = 0; gamm(12,5) = 0; ...
                gamm(12,6) = 0; gamm(12,7) = 5/2; ...
                gamm(12,8) = 16/315 - 160/19683*5/2; ...
                gamm(12,9) = 243/1540 + 35/2673*5/2; ...
                gamm(12,10) = 243/3080 + 7/2673*5/2; ...
                gamm(12,11) = 25/1386*(5 + sqrt15) - 3500/216513*(31 + 8*sqrt15)*5/2; ...
                gamm(12,12) = 25/1386*(5 - sqrt15) - 3500/216513*(31 - 8*sqrt15)*5/2;
            gamm(13,1) = 9/280; gamm(13,2) = 0; gamm(13,3) = 0; ...
                gamm(13,4) = 0; gamm(13,5) = 0; gamm(13,6) = 0; ...
                gamm(13,7) = 0; gamm(13,8) = 16/315; ...
                gamm(13,9) = 243/1540; gamm(13,10) = 243/3080; ...
                gamm(13,11) = 25/1386*(5+sqrt15); ...
                gamm(13,12) = 25/1386*(5-sqrt15); gamm(13,13) = 0;
            
            cdot = zeros(1,13);
            cdot(1) = 9/280;
            cdot(2) = 0;
            cdot(3) = 0;
            cdot(4) = 0;
            cdot(5) = 0;
            cdot(6) = 0;
            cdot(7) = 0;
            cdot(8) = 32/315;
            cdot(9) = 729/3080;
            cdot(10) = 729/3080;
            cdot(11) = 125/693;
            cdot(12) = 125/693;
            cdot(13) = 9/280;
            
            % RKN8(9)-13-\dot{x}
            f = zeros(13,3,length(bodies));
            rk = zeros(length(bodies), 3);
            for k = 1:13
                newsys = copy2(bodies,rk);
                for i = 1:length(bodies)
                    f(k,:,i) = getForce(newsys(i), newsys);
                    rk(i,:) = bodies(i).vel*alpha(k)*delt + gamm(k,1:k)*f(1:k,:,i)*delt^2;
                end
            end
                
            % get final values and update stats
            for i = 1:length(bodies)
                bodies(i).pos = bodies(i).pos + rk(i,:);
                bodies(i).vel = bodies(i).vel + cdot(1:13)*f(1:13,:,i)*delt;
            end
        end
    end
end
