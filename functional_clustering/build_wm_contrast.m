function [WMdiff, contrast_labels] = build_wm_contrast(WM)
%
% INPUT
%   WM : [regions x time x subjects x 8]
%
% OUTPUT
%   WMdiff : [regions x time x subjects x 4]
%            contrasts: 2bk - 0bk within category
%   contrast_labels : {'body','faces','places','tools'}

assert(ndims(WM) == 4, 'WM must be [regions x time x subjects x 8].');
assert(size(WM, 4) == 8, 'WM must have 8 stimulus slots.');

idx_0 = [1 2 3 4];
idx_2 = [5 6 7 8];

WMdiff = WM(:, :, :, idx_2) - WM(:, :, :, idx_0);
contrast_labels = {'body', 'faces', 'places', 'tools'};

end