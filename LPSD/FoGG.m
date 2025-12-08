function [Ga,Go,Gx,Gy] = FoGG(im,sigma,radius)
% Fast oriented gradients using derivative of Gaussian and atan2d for robust angle computation
GS = radius;
X = -GS:GS;
Y = -GS:GS;
[XX,YY] = meshgrid(X,Y);

G = (1 / (2 * pi * sigma^2)) * exp(-(XX.^2 + YY.^2) / (2 * sigma^2));
GxK = -XX / sigma^2 .* G;
GyK = -YY / sigma^2 .* G;

Gx = conv2(im,GxK,'same');
Gy = conv2(im,GyK,'same');

Ga = sqrt(Gx.^2 + Gy.^2);
% Use atan2d for correct quadrants and fewer conditionals
Go = atan2d(Gy, Gx);
Go(Go < 0) = Go(Go < 0) + 360; % map to [0,360)

