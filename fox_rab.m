function fox_rab(task)
% fox_rab solves a system of ODEs to determine if a fox can catch a rabbit 
% before it reaches it's burrow, within a specific area. Certain conditions
% about their running speeds, starting positions, and running paths are
% predetermined. fox_rab can take a single input parameter: 
% if the input is 1, then the constant speed model is used, 
% if the input is 2, then the diminishing speed model is used.

z0 = [-250, -550, 0, 0, 0, 0]; % Initial conditions of fox and rabbit positions and distances travelled
s_r0 = 13; s_f0 = 16; % Initial speeds of rabbit and fox respectively    
burrow = [-600 600];
mindist = 0.1;
S = [-200; -400];
tspan = [0 norm([z0(1) z0(2)] - [z0(4) z0(5)])/(s_f0 - s_r0)]; % Timespan to run ODE over

%ODE solver
option = odeset('MaxStep',0.01,'Events',@(t,z)events(z,burrow,mindist));
[~, pfoxrab, te, ze, ie] = ode45(@(t,z)fox_rab_ODE(z,s_f0,s_r0,task,S), tspan, z0, option);

xfox = pfoxrab(:,1);
yfox = pfoxrab(:,2);
xrab = pfoxrab(:,4);
yrab = pfoxrab(:,5);

% Plots of fox and rabbit paths, and warehouse for context
plot(xfox, yfox, xrab, yrab); hold on 
rectangle('Position',[-1000 -400 800 400]), text(-400,-200,'Warehouse')
legend('Fox', 'Rabbit')
xlim([-650 0]), ylim([-600 650])
xlabel('x coordinate (meters)'), ylabel('y coordinate (meters)')
hold off

% Display outputs with context

if ie == 1
    disp(['Rabbit reached the burrow before it was caught, in ', num2str(te(1)), ' seconds.' ...
        ' At this point the fox was at (', num2str(ze(1)), ', ', num2str(ze(2)), ').' ])
elseif ie == 2
    disp(['Rabbit was caught in ', num2str(te(1)), ' seconds at ' ...
        '(', num2str(ze(1)), ', ', num2str(ze(2)), ').' ])
end
disp(['The fox travelled ', num2str(ze(3)), ' meters in the ', num2str(te(1)), ' seconds.'])

    function dzdt = fox_rab_ODE(z,s_f0,s_r0, task, S)
        % Funtion sets the system of ODES for
        % position vector of and distance travelled
        % by fox and rabbit.
    
        dzdt = zeros(6,1);

        if task == 2
            % Diminishing speeds
            s_f = s_f0 * exp(-0.0002*z(3)); 
            s_r = s_r0 * exp(-0.0008*z(6));
        else
            % Constant speeds
            s_f = s_f0; s_r = s_r0; 
        end
    
        dist_rf = sqrt((z(4)-z(1))^2+(z(5)-z(2))^2); % Calc distance from rabbit to fox
        dist_Sf = sqrt((S(1)-z(1))^2 + (S(2)-z(2))^2); % Calc distance from fox to S
    
        if ~InSight(z) && z(1) < -200 % S blocks view of rabbit
            dzdt(1) = s_f/dist_Sf .* (S(1)-z(1));
            dzdt(2) = s_f/dist_Sf .* (S(2)-z(2));
        elseif ~InSight(z) % N block view of rabbit
            dzdt(1) = s_f .* 0;
            dzdt(2) = s_f .* 1;
        else % View of rabbit is not blocked
            dzdt(1) = s_f/dist_rf .* (z(4)-z(1));
            dzdt(2) = s_f/dist_rf .* (z(5)-z(2));
        end
        dzdt(3) = sqrt(dzdt(1)^2+dzdt(2)^2); % Distance travelled by fox at time t
    
        % ODE for rabbit's path
        dzdt(4) = -s_r/sqrt(2);
        dzdt(5) = s_r/sqrt(2);
        dzdt(6) = sqrt(dzdt(4)^2+dzdt(5)^2); % Distance travelled by rabbit at time t

    end
    
    function [value,isterminal,direction] = events(z, burrow, mindist)
    % Events function to stop ODE solver if rabbit reaches the burrow
    % or if fox catches the rabbit
    
    
    % Rabbit reaches the burrow x-coord
    value(1) = z(4) - burrow(1);
    isterminal(1) = 1;
    direction(1) = -1;


    % Fox catches the rabbit
    value(2) = sqrt((z(4)-z(1))^2+(z(5)-z(2))^2) - mindist;
    isterminal(2) = 1;
    direction(2) = -1;

    end

end