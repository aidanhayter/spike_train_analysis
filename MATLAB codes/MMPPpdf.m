%MMPPpdf plots the pdf of intervals (analytic form of the IIH) for the MMPP 
%process.  It takes as input the time axis length in seconds, and the number of 
%bins in it, as well as the four MMPP parameters.

trans1 = 101;
trans2 = 50;
lambda1 = 30;
lambda2 = 30;
tLength = 0.5;
nBins = 1000;
intLength = tLength/nBins;
a = trans1;
b = trans2;
c1 = lambda1;
c2 = lambda2;
t = [tLength/nBins:intLength:tLength];

A = a^2 + 2*a*c1 + 2*a*b -2*c2*a + c1^2 - 2*c1*b - 2*c1*c2 + b^2 +2*b*c2 + c2^2;
B = a + c1 + b +c2 + sqrt(A);
C = a + c1 + b + c2 - sqrt(A);

L1 = c1^2*b^2*sqrt(A) + 2*c1*b*c2*sqrt(A)*a + c2^3*a*c1 - c2^2*a^2*b;
L2 = -2*c2^2*a^2*c1 - c1^2*b^2*a + c1^3*b^2 + c2^3*a^2 - c2^2*a^3 - c1^2*b^3;
L3 = c1^2*b*c2*sqrt(A) - 2*c1*b*a^2*c2 + c1^3*b*c2 - c1^2*b*c2^2 - 2*c1*b^2*c2*a;
L4 = -c2*a*c1^2*b - c2^2*a*c1^2 - 2*c1^2*b^2*c2 - c1*b*c2^2*a + c2^2*a^2*sqrt(A);
L5 = c2^2*a*sqrt(A)*c1;

L6 = c1^2*b^2*sqrt(A) + c1^2*b*c2*sqrt(A) - c1^3*b^2 + c1^2*b^3 - c2^3*a^2 + c2^2*a^3 + 2*c2^2*a^2*c1;
L7 = 2*c1*b*c2*sqrt(A)*a - c1^3*b*c2 + 2*c1*b*a^2*c2 + c1*b*c2*a + c2^2*a*sqrt(A)*c1;
L8 = c1^2*b*c2^2 + c2^2*a^2*b + c2^2*a*c1^2 + c2*a*c1^2*b - c2^3*a*c1 + c1^2*b^2*a;
L9 = 2*c1^2*b^2*c2 + 2*c1*b^2*c2*a + c2^2*a^2*sqrt(A);

part1 = L1 + L2 + L3 + L4 + L5;
part2 = L6 + L7 + L8 + L9;
denom = (c1*b + c2*a)*C*sqrt(A)*B;

prob = -2*( part1*exp(-0.5*B*t)*(-0.5*B) + part2*exp(-0.5*C*t)*(-0.5*C) )/denom;
figure(2),
plot(t,prob);