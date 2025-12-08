function [kpts1,kpts2] = kpts_elimination(kpts1,kpts2,H,t)
% Eliminate keypoints that have no corresponding point within threshold under H and H^{-1}
pt1_H = H * [kpts1(:,1:2)'; ones(1,size(kpts1,1))];
pt1_H = (pt1_H(1:2,:) ./ pt1_H(3,:))';

pt2_H = (H \ [kpts2(:,1:2)'; ones(1,size(kpts2,1))]);
pt2_H = (pt2_H(1:2,:) ./ pt2_H(3,:))';

% Vectorized distance checks
if isempty(kpts1) || isempty(kpts2)
    kpts1 = []; kpts2 = [];
    return;
end
% Compute distances from each transformed kpts2 to all kpts1 and check threshold
dx2 = pt2_H(:,1) - kpts1(:,1)';
dy2 = pt2_H(:,2) - kpts1(:,2)';
D21min = sqrt(min(dx2.^2 + dy2.^2, [], 2));
idx2 = D21min < t;
kpts2 = kpts2(idx2, :);

if isempty(kpts2)
    kpts1 = [];
    return;
end
dx1 = pt1_H(:,1) - kpts2(:,1)';
dy1 = pt1_H(:,2) - kpts2(:,2)';
D12min = sqrt(min(dx1.^2 + dy1.^2, [], 2));
idx1 = D12min < t;
kpts1 = kpts1(idx1, :);
