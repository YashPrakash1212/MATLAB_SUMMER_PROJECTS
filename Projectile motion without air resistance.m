% -------------------------------------------------------------
% Projectile Motion Model (No Air Resistance)
% -------------------------------------------------------------
%
% In projectile motion, an object moves in two directions:
%
% 1. Horizontal motion (x-direction)
%    - There is no acceleration (ignoring air resistance)
%    - Velocity stays constant
%
% 2. Vertical motion (y-direction)
%    - Gravity causes a constant downward acceleration
%    - The object slows down while going up
%    - The object speeds up while falling down
%
% The motion is calculated using kinematic equations:
%
% final_velocity = initial_velocity + acceleration*time
%
% final_position = initial_position + initial_velocity*time
%                  + 0.5*acceleration*time^2


clc
clear


% -------------------------------------------------------------
% Initial Conditions
% -------------------------------------------------------------

speed = 50;             
% Initial speed of the projectile (m/s)
% This is the speed at the exact moment the object is launched.


angle_of_launch = 45;   
% Angle that the projectile is launched above the horizontal (degrees)


gravity = 9.8;          
% Acceleration caused by Earth's gravity (m/s^2)
% Gravity always acts downward, so it is negative in the equation.


initial_vertical_position = 0; 
% Starting height of projectile (meters)
% We assume it starts at ground level.


% -------------------------------------------------------------
% Breaking Velocity Into Components
% -------------------------------------------------------------

initial_vertical_velocity = speed*sind(angle_of_launch);

% The vertical velocity controls:
% - how high the projectile goes
% - how long it stays in the air
%
% We use sine because we are finding the vertical component.


initial_horizontal_velocity = speed*cosd(angle_of_launch);

% The horizontal velocity controls:
% - how far the projectile travels
%
% We use cosine because we are finding the horizontal component.


% -------------------------------------------------------------
% Finding Time of Flight
% -------------------------------------------------------------

time_to_max_height = initial_vertical_velocity/gravity;

% At maximum height:
% vertical velocity becomes zero.
%
% Using:
% final_velocity = initial_velocity + acceleration*time
%
% 0 = initial_vertical_velocity - gravity*time
%
% Solving gives:
% time = initial_vertical_velocity/gravity


total_time = 2*time_to_max_height;

% The projectile takes the same amount of time going down
% as it takes going up because it lands at the same height.


% -------------------------------------------------------------
% Finding Maximum Height
% -------------------------------------------------------------

maximum_height = initial_vertical_position ...
    + (initial_vertical_velocity*time_to_max_height) ...
    - (0.5*gravity*time_to_max_height^2);

% This uses:
%
% y = yi + v0y*t - 0.5*g*t^2
%
% We use the time at maximum height because that is when
% the projectile reaches its highest point.


% -------------------------------------------------------------
% Finding Range
% -------------------------------------------------------------

range = initial_horizontal_velocity*total_time;

% Horizontal motion has constant velocity:
%
% distance = velocity*time
%
% Range tells us how far horizontally the projectile travels.


% -------------------------------------------------------------
% Display Results
% -------------------------------------------------------------

fprintf("Maximum Height: %.2f meters\n", maximum_height)

fprintf("Time of Flight: %.2f seconds\n", total_time)

fprintf("Range: %.2f meters\n", range)



% -------------------------------------------------------------
% Plotting Height vs Time
% -------------------------------------------------------------

time = linspace(0,total_time,100);

% Creates 100 time points between launch and landing.


height = initial_vertical_position ...
    + (initial_vertical_velocity*time) ...
    - (0.5*gravity*time.^2);

% Calculates the height of the projectile at each time.


figure

plot(time,height,'LineWidth',2)

grid on

xlabel('Time (seconds)')

ylabel('Height (meters)')

title('Projectile Motion: Height vs Time')