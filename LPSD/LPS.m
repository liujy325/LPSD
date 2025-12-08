function LPP = LPS(im,r_num,t_num,theta)
theta = mod(theta,360);
[row,col] = size(im);
LPP = zeros(t_num,r_num);

Rmax = min(row,col)/2-1;
Rmin = 7;
if r_num>1
    dR = (Rmax/Rmin)^(1/(r_num-1));
else
    dR = 1;
end
dt = 2*pi/t_num;
RI = Rmin * dR.^(0:r_num-1);
TI = (0:t_num-1) * dt;

k = mod(round(theta/360*t_num), t_num);

% Center coordinates
y0 = round(row/2); x0 = round(col/2);

% Precompute cos/sin tables
ct = cos(TI);
st = sin(TI);

% Generate all sample coordinates using outer sums
[Rgrid,Tgrid] = ndgrid(RI, 1:t_num); % Rgrid: r_num x t_num, Tgrid indexes angles
xs = round(x0 + Rgrid .* ct(Tgrid));
ys = round(y0 + Rgrid .* st(Tgrid));

xs = min(max(xs,1), col);
ys = min(max(ys,1), row);

% Linear indexing and gather values
linIdx = sub2ind([row,col], ys(:), xs(:));
vals = im(linIdx);
vals = reshape(vals, r_num, t_num).'; % -> t_num x r_num

if k>0
    LPP(1:t_num-k,:) = vals(k+1:end,:);
    LPP(t_num-k+1:end,:) = vals(1:k,:);
else
    LPP = vals;
end
end