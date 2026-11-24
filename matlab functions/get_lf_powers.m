
function [lf_powers, entire_lf_within_coi] = get_lf_powers(bp_uniform, fs)


bp_uniform = bp_uniform(:);

% run cwt
%cwt(bp_uniform, fs, "VoicesPerOctave",16); % TEST: plot scalogram
[wt, frequencies, coi] = cwt(bp_uniform, fs, "VoicesPerOctave",16);

% reorder coefs and frequencies
[sorted_frequencies, order] = sort(frequencies);

sorted_wt = wt(order, :);

% calculate power matrix
power_matrix = abs(sorted_wt).^2;

% create coi and lf combined mask
coi_mask = sorted_frequencies(:) >= coi(:)';
lf_mask = sorted_frequencies >= 0.04 & sorted_frequencies <= 0.15;
lf_coi_mask = coi_mask(lf_mask, :);

entire_lf_within_coi = all(lf_coi_mask, 1);

% calculate LF power
N = size(power_matrix, 2);
lf_powers = zeros(1, N);
for i = 1:N
    timeframe_powers = power_matrix(:, i);
    lf_power = trapz(sorted_frequencies(lf_mask), timeframe_powers(lf_mask));
    lf_powers(i) = lf_power;
end

end
