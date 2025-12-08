function [kpts,im] = detector(im,N)
% ETM and GridFAST with faster orientation estimation
[Ga1,~] = FoGG(im,2.5,3);
[Ga2,~] = FoGG(Ga1,2.5,3);
for k = 1:5
    Ga1 = max(Ga1 - Ga2, 0);
    [Ga2,~] = FoGG(Ga1,2.5,3);
end
im = Ga1;

kpts = gridFAST(im,N,30,0.001,0.001);

end

function res = gridFAST(im,N,gridnum,MinContrast,MinQuality)
% N: k_max
% gridnum: im is divided into gridnum*gridnum regions
key = detectFASTFeatures(im,"MinContrast",MinContrast,"MinQuality",MinQuality);
location = key.Location;
if size(location,1) <= N
    res = [key.Location ones(size(key.Location,1),1) zeros(size(key.Location,1),1)];
    return
end

grid_size = double(size(im))/gridnum;
grid = floor(location(:,1)/grid_size(2)) + gridnum*floor(location(:,2)/grid_size(1))  + 1;   % grid ID
kpts = [key.Location ones(size(key.Location,1),1) zeros(size(key.Location,1),1) key.Metric grid];   % kpts is N*6 matrix: x,y,scale,angle,metric,grid
kpts = sortrows(kpts,5,"descend");
kpts_subset = cell(gridnum,gridnum);

for i = 0:gridnum-1
    for j = 0:gridnum-1
        kpts_subset(j+1,i+1) = {kpts(kpts(:,6)==i+j*gridnum+1,:)};
    end
end

if size(location,1) < N
    N = size(location,1);
end
kpts_C = zeros(N,6);
cnt = 0;
while cnt < N
    for i = 1:gridnum^2
        if cnt >= N
            break;
        end
        if ~isempty(kpts_subset{i})
            cnt = cnt + 1;
            kpts_tmp = kpts_subset{i};
            kpts_C(cnt,:) = kpts_tmp(1,:);
            kpts_tmp(1,:) = [];
            kpts_subset{i} = kpts_tmp;
        end
    end
end

res = kpts_C(:,1:4);
end