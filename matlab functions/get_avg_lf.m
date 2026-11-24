% TO DO:
% make it so that you can customize the time window of analysis - done
% connect the real bp python data to the matlab file - done
% figure out if they only put the SBP values in here - done
% verify if this whole thing is even done correctly????? - done
% make plots!!!

function [avg_lf_power] = get_avg_lf(lf_powers, sample_seconds, fs, entire_lf_within_coi, start_idx, end_idx) % temporarily remove input signal

%error handling
if sample_seconds * fs > size(lf_powers, 2)
    error("Analysis window is longer than the recording.");
end

% check if the middle 40 seconds is within the COI
analysis_sample = lf_powers(start_idx:end_idx);

valid_mask = entire_lf_within_coi(start_idx : end_idx);

valid_analysis_sample = analysis_sample(valid_mask);

if length(valid_analysis_sample) < (fs * sample_seconds) - (fs * 5)
    disp("Warning: Possibly poor analysis results")
    
end

disp(length(valid_analysis_sample) + " valid data points out of " + length(analysis_sample) + " data points from analysis window" )
% get average over analysis window
avg_lf_power = mean(valid_analysis_sample);

end



