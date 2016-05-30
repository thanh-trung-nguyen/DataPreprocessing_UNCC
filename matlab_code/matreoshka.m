function [x,y] = matreoshka



dt = pi/200;
t = -pi/2:dt:pi/2;

Nt = length(t);
x = zeros(1,Nt); y = x;


Angle0 = acos(2.5/3.9);
idx1 = find(t < -Angle0);
x(idx1) = linspace(0,0.03,length(idx1));
y(idx1) = -0.037;

idx2 = find(t >= -Angle0);
idx3 = find(t <= Angle0);
idx2 = idx2(1):idx3(end);

x(idx2) = 0.039*cos(t(idx2));
y(idx2) = 0.039*sin(t(idx2));

Angle1 = acos(2.5/sqrt(2.5^2 + 6^2));
idx3 = find(t > Angle0);
idx4 = find(t < Angle1);
idx3 = idx3(1): idx4(end);

x(idx3) = 0.025; 
y(idx3) = linspace(0.03,0.06,length(idx3));

idx4 = find(t >=Angle1);
t1 = linspace(acos(2.5/3.6),pi/2,length(idx4));

x(idx4) = 0.036*cos(t1);
y(idx4) = 0.034 + 0.036*sin(t1);

x = [x, -x(end:-1:1)];
y = [y, y(end:-1:1)];



