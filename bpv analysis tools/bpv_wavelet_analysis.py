from process_data import timestamp_extraction, low_pass_filter, extract_peaks, remove_artifacts, interpolate_missed_beats, interpolate_and_resample, extract_subject_from_filename, manual_peak_detection, hhmmss_to_seconds
from bpv_metrics import get_lf_powers, get_avg_lf
import numpy as np
import pandas as pd
import os
import math
import matplotlib.pyplot as plt
from pathlib import Path


# ############################################# USER Configuration ##################################################
SAMPLING_RATE = 1000
SAMPLE_SECONDS = 40


RESAMPLED_RATE = 2
PROMINENCE = 4
WINDOW_SIZE = 20
UPPER_LIMIT = 1.2

############################################################################################################################

def file_processing():

    string_bp_file = input("Enter relative filepath of bp data file: ")
    bp_file = Path(string_bp_file)

    string_output_folder = input("enter output folder relative filepath: ")

    output_folder = Path(string_output_folder)

    if bp_file.exists() and bp_file.is_file():
        t, peaks = analyze_bpv(bp_file)

    # get timestamps
    start_baseline = input("Pre stim window start timestamp (hh:mm:ss): ") 
    end_baseline = input("Pre stim window end timestamp (hh:mm:ss): ")

    start_during = input("During stim start timestamp (hh:mm:ss): ")
    end_during = input("During stim end timestamp (hh:mm:ss): ")

    start_post = input("Post stim start timestamp (hh:mm:ss): ")
    end_post = input("Post stim end timestamp (hh:mm:ss): ")

    # convert timestamps to samples
    start_baseline_samples = (hhmmss_to_seconds(start_baseline) + t[0]) * RESAMPLED_RATE# start total signal: 3480
    end_baseline_samples = hhmmss_to_seconds(end_baseline) * RESAMPLED_RATE # end baseline: 3660

    start_during_samples = hhmmss_to_seconds(start_during) * RESAMPLED_RATE
    end_during_samples = hhmmss_to_seconds(end_during) * RESAMPLED_RATE # end during: 3928

    start_post_samples = hhmmss_to_seconds(start_post) * RESAMPLED_RATE
    end_post_samples = len(peaks)
    

    # normalize samples in terms of starting peak
    end_baseline_samples = end_baseline_samples - start_baseline_samples

    start_during_samples = start_during_samples - start_baseline_samples
    end_during_samples = end_during_samples - start_baseline_samples

    start_post_samples = start_post_samples - start_baseline_samples
    start_baseline_samples = 0

    # get start and end idxs for middle window
    n = RESAMPLED_RATE * SAMPLE_SECONDS

    baseline_start_idx = math.ceil((end_baseline_samples + start_baseline_samples)/2) - math.floor(n/2)
    baseline_end_idx = baseline_start_idx + n - 1
        
    during_start_idx = math.ceil((end_during_samples + start_during_samples)/2) - math.floor(n/2)
    during_end_idx = during_start_idx + n - 1

    post_start_idx = math.ceil((end_post_samples + start_post_samples)/2) - math.floor(n/2)
    post_end_idx = post_start_idx + n - 1

    # # Calculate LF powers
    lf_powers, entire_lf_within_coi = get_lf_powers(peaks, RESAMPLED_RATE)

    # Calculate average lf powers within windows
    baseline_avg_lf = get_avg_lf(lf_powers, SAMPLE_SECONDS, RESAMPLED_RATE, entire_lf_within_coi, baseline_start_idx, baseline_end_idx)
    during_avg_lf = get_avg_lf(lf_powers, SAMPLE_SECONDS, RESAMPLED_RATE, entire_lf_within_coi, during_start_idx, during_end_idx)
    post_avg_lf = get_avg_lf(lf_powers, SAMPLE_SECONDS, RESAMPLED_RATE, entire_lf_within_coi, post_start_idx, post_end_idx)

    bpv_full_df = pd.DataFrame({"bpv parameters" : ["Average Wavelet LF Power"], "Before" : [baseline_avg_lf], "During" : [during_avg_lf], "After" : [post_avg_lf], "Delta" : [post_avg_lf - baseline_avg_lf]})

    subject = extract_subject_from_filename(str(bp_file))
    #folder_path = os.path.join(output_folder)
    os.makedirs(output_folder, exist_ok=True)
    output_filename = output_folder / f"{subject}-WAVELET_bpv.csv"
    bpv_full_df.to_csv(output_filename, index=False)

def analyze_bpv(input_file):
    bp = np.loadtxt(input_file)

    bp = timestamp_extraction(bp, SAMPLING_RATE)

    ts = 1/SAMPLING_RATE
    t = np.arange(0, len(bp)/SAMPLING_RATE, ts)

    bp = low_pass_filter(bp, SAMPLING_RATE)

    peak_idxs = extract_peaks(bp, t, PROMINENCE)

    peak_vals, peak_idxs = remove_artifacts(bp, t, peak_idxs, WINDOW_SIZE)
   
    manual_peak_detection(bp, t, peak_idxs, bp[peak_idxs])

    peak_vals, peak_times = interpolate_missed_beats(bp, t, bp[peak_idxs], peak_idxs, WINDOW_SIZE, UPPER_LIMIT)

    peak_vals_uniform, t_uniform = interpolate_and_resample(peak_vals, peak_times, RESAMPLED_RATE)
    
    # plt.plot(t, bp, color="orange")
    # plt.scatter(t_uniform, peak_vals_uniform, color="blue")
    # plt.title("Resampled BP with plotted Systolic Peaks")
    # plt.show()

    return t_uniform, peak_vals_uniform

####################################################################################################3
######################################## RUN HERE ###################################################
file_processing()


