function [des,kpts] = descriptor_LPS_sen(im,kpts,patch_size,nblock,nbin,r_num,t_num,H)
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

H(abs(H)<1e-4) = 0;
[Q,~] = qr(H(1:2,1:2));
ang = atan2d(Q(2,1),Q(1,1));
if H(1,1) > 0
    ang = ang + 180;
end
angidx = round(mod(ang,180)/(180/nbin));
% Adjust bins globally once
if angidx~=0
    QOM(QOM>0) = QOM(QOM>0) + nbin - angidx;
    QOM(QOM>nbin) = QOM(QOM>nbin) - nbin;
end

for k = 1:n
    x = X(k); y = Y(k);
    [qom_LPS,isout] = LPS_H(QOM,x,y,r,r_num,t_num,H);
    if isout
        continue
    else
        histo = LPSD(qom_LPS,nblock,nbin);
    end
    nrm = norm(histo);
    if nrm>0
        des(:,k) = histo / nrm;
    end
end

idx = any(des);
des = des(:,idx);
kpts = kpts(:,idx)';

end