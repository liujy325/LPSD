function [LPP,isout] = LPS_H(im,x,y,r,r_num,t_num,H)
% Homography-aware Log-Polar Sampling
isout = 0;
x0 = x + r;
y0 = y + r;
[row,col] = size(im);
LPP = zeros(t_num,r_num);

Rmax = r;
Rmin = 7;
if r_num>1
    dR = (Rmax/Rmin)^(1/(r_num-1));
else
    dR = 1;
end
dt = 2*pi/t_num;
RI = Rmin * dR.^(0:r_num-1);
TI = (0:t_num-1) * dt;

% Only affine part of H is used per original code
H2 = [H(1,1),H(1,2),0; H(2,1),H(2,2),0; 0,0,1];

% Precompute cos/sin
ct = cos(TI);
st = sin(TI);
[Rgrid,Tgrid] = ndgrid(RI, 1:t_num); % r_num x t_num
dx = Rgrid .* ct(Tgrid);
dy = Rgrid .* st(Tgrid);

% Apply transform to all offsets
onesGrid = ones(size(dx));
Xw = H2(1,1).*dx + H2(1,2).*dy; % no translation term
Yw = H2(2,1).*dx + H2(2,2).*dy;
W  = onesGrid; % third row is 1
Xn = Xw./W; Yn = Yw./W;

xs = round(x0 + Xn);
ys = round(y0 + Yn);

% Bounds check; if any sample out of bounds, mark isout and return
if any(xs(:)<=0 | xs(:)>col | ys(:)<=0 | ys(:)>row)
    isout = 1;
    return;
end

linIdx = sub2ind([row,col], ys(:), xs(:));
vals = im(linIdx);
LPP = reshape(vals, r_num, t_num).'; % t_num x r_num
end
