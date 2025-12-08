function [des,kpts] = descriptor_LPS(im,kpts,patch_size,nblock,nbin,r_num,t_num)
pad = floor(patch_size/2);
im = padarray(im,[pad,pad]);
[Ga, Go] = FoGG(im,2.5,3);
QOM = floor(mod(Go,180) ./ (180/nbin)) + 1;
QFoGG = floor(mod(Go,360) ./ (180/nbin)) + 1;

GaMin = min(Ga,[],"all"); GaMax = max(Ga,[],"all");
if GaMax>GaMin
    GaN = (Ga-GaMin)/(GaMax-GaMin);
else
    GaN = zeros(size(Ga), 'like', Ga);
end
Tmin = 0.05;
QOM(GaN<Tmin) = 0;
QFoGG(GaN<Tmin) = 0;

n = size(kpts,2);
des = zeros(nblock*nblock*nbin,n);
X = round(kpts(1,:));
Y = round(kpts(2,:));
r = floor(patch_size/2);

for k = 1:n
    x = X(k); y = Y(k);
    x = x + r; y = y + r;
    qom = QOM(y-r:y+r-1, x-r:x+r-1);
    ip = dominant_orientation(QOM,QFoGG,x,y,round(r/2),nbin);
    ip_norm = ip;
    ip_norm(ip_norm>nbin) = ip_norm(ip_norm>nbin)-nbin;
    if any(qom(:))
        qom = qom - (ip_norm - nbin).*(qom>0);
        over = qom>nbin; qom(over) = qom(over) - nbin;
    end
    qom_LPS = LPS(qom,r_num,t_num,(ip-1)/nbin*180);
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

function ip = dominant_orientation(QOM,QFoGG,x,y,r,nbin)
    qom = QOM(y-r:y+r-1, x-r:x+r-1);
    qfogg = QFoGG(y-r:y+r-1, x-r:x+r-1);
    h1 = histcounts(qom(:),nbin+1);
    h2 = histcounts(qfogg(:),2*nbin+1);
    [~,ip] = max(h1(2:end));
    ischg = sum(h2(2:nbin+1)) < sum(h2(nbin+2:2*nbin+1));
    if ischg
        ip = ip + nbin;
    end
end