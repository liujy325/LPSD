function [des,kpts] = descriptor_LPS_ref(im,kpts,patch_size,nblock,nbin,r_num,t_num)
pad = floor(patch_size/2);
im = padarray(im,[pad,pad]);
[Ga, Go] = FoGG(im,2.5,3);
QOM = floor(mod(Go,180)./(180/nbin)) + 1;
GaMin = min(Ga,[],"all"); GaMax = max(Ga,[],"all");
if GaMax>GaMin
    GaN = (Ga-GaMin)/(GaMax-GaMin);
else
    GaN = zeros(size(Ga), 'like', Ga);
end
Tmin = 0.05;
QOM(GaN<=Tmin) = 0;

n = size(kpts,2);
des = zeros(nblock*nblock*nbin,n);
X = round(kpts(1,:));
Y = round(kpts(2,:));
r = floor(patch_size/2);

for k = 1:n
    x = X(k); y = Y(k);
    x = x + r; y = y + r;
    qom = QOM(y-r:y+r-1, x-r:x+r-1);
    qom_LPS = LPS_R(qom,r_num,t_num,0);
    histo = LPSD(qom_LPS,nblock,nbin);
    nrm = norm(histo);
    if nrm>0
        des(:,k) = histo / nrm;
    end
end

idx = any(des);
des = des(:,idx);
kpts = kpts(:,idx)';

end