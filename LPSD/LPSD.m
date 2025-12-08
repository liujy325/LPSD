function histo = LPSD(pim, ns, nbin)
% Log-Polar Sector Descriptor histogram per ns x ns grid, vectorized
[h, w] = size(pim);
grid_h = floor(h / ns); % rows per block
grid_w = floor(w / ns); % cols per block
h = grid_h * ns; w = grid_w * ns;
pim = pim(1:h, 1:w);

% Block indices for each pixel
[ii,jj] = ndgrid(1:h, 1:w);
bi = ceil(ii / grid_h); % 1..ns
bj = ceil(jj / grid_w); % 1..ns
blockId = (bi-1)*ns + bj; % 1..ns*ns

% Only positive bins contribute
mask = pim > 0;
vals = pim(mask);
bids = blockId(mask);

% Accumulate counts per (block, bin)
% Linearize (block,bin) -> idx = (bin-1)*ns*ns + block
idx = (vals-1)*(ns*ns) + bids; % vals in 1..nbin
counts = accumarray(idx, 1, [nbin*ns*ns, 1]);
histo = reshape(counts, [ns*ns, nbin]).'; % nbin x (ns*ns)

% Normalize per block area
histo = histo / (grid_h * grid_w);
histo = histo(:);
end
