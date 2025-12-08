function [solution,rmse,cor1,cor2]= LPSD_match(des1,loc1,des2,loc2,change_form,error_t,isunique)
if (size(loc1,1) < 3 || size(loc2,1) < 3)
    solution = [];
    rmse = inf;
    cor1 = [];
    cor2 = [];
    return
end
%% NNDR Match
distRatio = 0.999;
Nmax=max(300,round(0.05*size(loc1,1)));
[cor11_low,cor22_low,cor11_up,cor22_up]= NNDR(des1,loc1,des2,loc2,distRatio,Nmax,isunique);

%% Remove Duplicate Points
uni1=[cor11_low(:,[1 2]),cor22_low(:,[1 2])];
[~,i,~]=unique(uni1,'rows','first');
cor11_low=cor11_low(sort(i)',1:2);
cor22_low=cor22_low(sort(i)',1:2);

if (size(cor11_low,1) < 3 || size(cor22_low,1) < 3)
    solution = [];
    rmse = inf;
    cor1 = [];
    cor2 = [];
    return
end

%% FSC
[solution,rmse,cor1,cor2]=FSC2(cor11_low,cor22_low,cor11_up,cor22_up,change_form,error_t);

end


function [cor11_low,cor22_low,cor11_up,cor22_up] = NNDR(des1,loc1,des2,loc2,distRatio,N,isunique)
[idx, dists] = knnsearch(des2, des1, 'K', 2);

tratios = dists(:,1) ./ dists(:,2);
tindxs  = idx(:,1);

[~, indx1] = sort(tratios);
tmatch1 = zeros(1, size(des1,1));
if length(indx1) < N
    tmatch1(indx1(1:end)) = tindxs(indx1(1:end));
else
    tmatch1(indx1(1:N)) = tindxs(indx1(1:N));
end

indx2 = find(tratios < distRatio);
tmatch2 = zeros(1, size(des1,1));
tmatch2(indx2) = tindxs(indx2);

[~, point1_low, point2_low] = find(tmatch1);
cor11_low = loc1(point1_low, 1:2);
cor22_low = loc2(point2_low, 1:2);

[~, i, ~] = unique(cor11_low, 'rows', 'first');
cor11_low = cor11_low(sort(i)', :);
cor22_low = cor22_low(sort(i)', :);
[~, i, ~] = unique(cor22_low, 'rows', 'first');
cor11_low = cor11_low(sort(i)', :);
cor22_low = cor22_low(sort(i)', :);

[~, point1_up, point2_up] = find(tmatch2);
cor11_up = loc1(point1_up, 1:2);
cor22_up = loc2(point2_up, 1:2);

if isunique
    [~, i, ~] = unique(cor11_up, 'rows', 'first');
    cor11_up = cor11_up(sort(i)', :);
    cor22_up = cor22_up(sort(i)', :);
    [~, i, ~] = unique(cor22_up, 'rows', 'first');
    cor11_up = cor11_up(sort(i)', :);
    cor22_up = cor22_up(sort(i)', :);
end
end
