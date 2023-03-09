a = -16;
b = -5;
DT = 0.0075*390;

tt = 0:8600:(365*86400);

figure
plot(tt/86400/365, a*cos(2*pi*tt/86400/365) + b + DT)
grid()

tstar1 = acos( (b + DT)/a)/2/pi;
tstar2 = acos(-(b + DT)/a)/2/pi;
xline(tstar1, 'r')
xline(tstar2, 'b')

