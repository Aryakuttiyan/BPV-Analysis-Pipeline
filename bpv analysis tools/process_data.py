import pandas as pd
from scipy.signal import find_peaks, butter,filtfilt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from scipy.interpolate import CubicSpline, PchipInterpolator
import os
import re

def low_pass_filter(bp, sampling_rate):
    
    cutoff = 10  
    order = 2 

    b, a = butter(order, cutoff,  btype='low', analog=False, fs=sampling_rate)
    filtered_signal = filtfilt(b, a, bp)

    t = np.arange(len(bp)) / sampling_rate

    return filtered_signal

def hhmmss_to_seconds(timestamp):
    return int(timestamp[:2]) * 3600 + int(timestamp[3:5]) * 60 + int(timestamp[6:])

def timestamp_extraction(data, sampling_rate):
    start = input("what is the starting timestamp? Please use hh:mm:ss format.")
    end = input("what is the end timestamp? Please use hh:mm:ss format.")
    start_seconds = hhmmss_to_seconds(start)
    end_seconds = hhmmss_to_seconds(end)

    if start_seconds == 0:
        start_idx = 0
        end_idx = end_seconds * sampling_rate
    else:
        start_idx = start_seconds * sampling_rate - 1
        end_idx = end_seconds * sampling_rate - 1
    data = data[start_idx : end_idx]

    return data

def extract_peaks(bp, t, prominence):

    peak_idxs, _ = find_peaks(bp, prominence=prominence)

    #Uncomment to see full BP waveform with SYS peaks
    # plt.plot(t, bp, color="orange")
    # plt.scatter(t[peak_idxs], bp[peak_idxs], color="blue")
    # plt.title("BP with plotted Systolic Peaks")
    # plt.show()
    return list(peak_idxs)

def remove_artifacts(bp, t, peak_idxs, window_size):
    
    peak_vals = drop_by_interval(bp, t, peak_idxs, window_size)
   
    peak_vals, peak_idxs = drop_nans(peak_vals, peak_idxs)

    peak_vals = drop_by_amplitude(bp, peak_idxs, window_size)
    peak_vals, peak_idxs = drop_nans(peak_vals, peak_idxs)
    return peak_vals, peak_idxs

def drop_nans(peak_vals, peak_idxs):

    nan_indices = np.where(np.isnan(peak_vals))[0]

    peak_vals = np.delete(peak_vals, nan_indices)
    peak_idxs = [
            x for i, x in enumerate(peak_idxs)
            if i not in list(nan_indices)
    ]

    return peak_vals, peak_idxs

def calculate_mads(x):
    median = np.median(x)
    mads = np.median(np.abs(x - median))
    return mads

def drop_by_interval(bp, t, peak_idxs, window_size):
    
    peak_vals = bp[peak_idxs]
    peak_times = t[peak_idxs]

    intervals = np.diff(peak_times)
    windows = pd.Series(intervals).rolling(window_size, center=True, min_periods=1)
    
    means = windows.apply(lambda x : np.mean(x))
    stds = windows.apply(lambda x : np.std(x))

    z_scores = np.where(stds == 0, np.nan, (intervals - means) / stds)

    for i in range(len(intervals) - 1):
        if abs(z_scores[i]) > 3 and abs(z_scores[i + 1]) > 3:
            peak_vals[i + 1] = np.nan

    return peak_vals

def interpolate_missed_beats(bp, t, peak_vals, peak_idxs, window_size, upper_limit):
    peak_times = t[peak_idxs]
    # experiment here: different interpolation methods (cublic spline, other stuff idk), base interpolationg on various number of points
    intervals = np.diff(peak_times)
    medians = pd.Series(intervals).rolling(window_size, center=True, min_periods=1).median().to_numpy() # change this

    interpolated_xs = np.array([])
    interpolated_ys = np.array([])
  
    peak_times_interpolated = peak_times
    peak_vals_interpolated = peak_vals
  
    for i in range(len(intervals)):
        
        if intervals[i] > medians[i] * upper_limit: # change this
            
            midpoint = (peak_times[i] + peak_times[i + 1])/2
            
            if i >= 5 and i + 7 <= len(peak_times):
                pchip = PchipInterpolator(peak_times[i - 5: i + 7], peak_vals[i - 5 : i + 7])
                y_interpolated = float(pchip(midpoint))
                
            else:
                pchip = PchipInterpolator([peak_times[i], peak_times[i + 1]], [peak_vals[i], peak_vals[i + 1]])
                y_interpolated = float(pchip(midpoint))
                
            interpolated_xs = np.append(interpolated_xs, midpoint)
            interpolated_ys = np.append(interpolated_ys, y_interpolated)
            time_idx = np.searchsorted(peak_times_interpolated, midpoint)
            peak_times_interpolated = np.insert(peak_times_interpolated, time_idx, midpoint)
            peak_vals_interpolated = np.insert(peak_vals_interpolated, time_idx, y_interpolated)

    # plt.plot(t, bp, color="orange")
    # plt.scatter(peak_times_interpolated, peak_vals_interpolated, color="blue")
    # plt.scatter(interpolated_xs, interpolated_ys, color="green", s=25)
    # plt.title("BP with plotted Systolic and Diastolic Peaks outliers filtered out and interpolated missing values")
    # plt.show()

    return peak_vals_interpolated, peak_times_interpolated
            
def drop_by_amplitude(bp, peak_idxs, window_size):
    peak_vals = bp[peak_idxs]

    windows = pd.Series(peak_vals).rolling(window_size, center=True, min_periods=1)
    
    means = windows.apply(lambda x : np.mean(x))
    stds = windows.apply(lambda x : np.std(x))

    z_scores = np.where(stds == 0, np.nan, (peak_vals - means) / stds)

    for i in range(len(peak_vals)):
        if abs(z_scores[i]) > 3:
            peak_vals[i] = np.nan

    return peak_vals

def interpolate_and_resample(peak_vals, peak_times, fs_new):
    t_uniform = np.arange(
    peak_times[0],
    peak_times[-1] + 1/fs_new,
    1/fs_new
    )

    cs = CubicSpline(peak_times, peak_vals)
    bp_uniform = cs(t_uniform)

    return bp_uniform, t_uniform

def extract_subject_from_filename(input_file):
    filename = os.path.basename(input_file)
    filename = os.path.splitext(filename)[0]
    return filename

def clean_filename(filename):
    filename = re.sub(r'[<>:"/\\|?*]', '', filename)
    filename = filename.strip(' .')
    return filename

def manual_peak_detection(bp, t, peak_idxs, peak_vals):

    
    user_deleted_peaks = []
    user_added_peaks = []
    fig, ax = plt.subplots()
    ax.plot(t, bp, color="orange")
    #plt.scatter(t[peak_idxs], bp[peak_idxs], color="green", s=60)
    points = ax.scatter(t[peak_idxs], peak_vals, color=["blue" for i in range(len(peak_vals))], picker=True)
    #plt.scatter(t[trough_idxs], trough_vals, color="green", s=10)
    ax.set_title("Click to delete false peaks")

    button_ax = fig.add_axes((0.7, 0.005, 0.1, 0.075))

    button = Button(button_ax, "Save")
    button.on_clicked(lambda x: save_deleted_changes(x, peak_idxs, t, user_deleted_peaks))
    fig.canvas.mpl_connect('pick_event', lambda x: delete_false_peaks(x, fig, points, user_deleted_peaks))
    plt.show()

    fig, ax = plt.subplots()
    ax.plot(t, bp, color="orange", picker=1)
    #plt.scatter(t[peak_idxs], bp[peak_idxs], color="green", s=60)
    points = ax.scatter(t[peak_idxs], bp[peak_idxs], color="blue")
    #plt.scatter(t[trough_idxs], trough_vals, color="green", s=10)
    ax.set_title("Add missed peaks")

    button_ax = fig.add_axes((0.7, 0.005, 0.1, 0.075))

    button = Button(button_ax, "Save")
    button.on_clicked(lambda x: save_added_changes(x, peak_idxs, t, user_added_peaks))
    fig.canvas.mpl_connect('pick_event', lambda x: add_missed_peak(x, fig, ax, user_added_peaks))
    plt.show()

def add_missed_peak(event, fig, ax, user_added_peaks):
    if event.mouseevent.inaxes != ax:
        return

    x_clicked = event.mouseevent.xdata
    y_clicked = event.mouseevent.ydata

    # Put a point at the clicked location
    ax.scatter(x_clicked, y_clicked, color="blue")
    user_added_peaks.append((x_clicked, y_clicked))
   
    fig.canvas.draw_idle()

def delete_false_peaks(event, fig, points, user_deleted_peaks):
    idx = event.ind[0]
    colors = points.get_facecolors()
    colors[idx] = [1, 0, 0, 1]
    points.set_facecolors(colors)
    fig.canvas.draw_idle()
    user_deleted_peaks.append(idx)

def save_deleted_changes(event, peak_idxs, t, user_deleted_peaks):
    if user_deleted_peaks != []:
        peak_idxs[:] = [
            x for i, x in enumerate(peak_idxs)
            if i not in user_deleted_peaks
        ]
    plt.close()

def save_added_changes(event, peak_idxs, t, user_added_peaks):
    if user_added_peaks != []:
        for coords in user_added_peaks:
            time_idx = np.searchsorted(t, coords[0])
            peak_idxs.append(time_idx)

    peak_idxs.sort()
    plt.close()



# bp = np.loadtxt(r"data\raw bp\VNS01_baseline_mapping_left_VNS_40_10_19_2026_bp.txt")
# bp = timestamp_extraction(bp, 1000)

# ts = 1/1000
# t = np.arange(0, len(bp)/1000, ts)

# bp = low_pass_filter(bp, 1000)

# peak_idxs = extract_peaks(bp, t, 4)

# peak_vals, peak_idxs = remove_artifacts2(bp, t, peak_idxs, 20)
# #trough_vals, trough_idxs = remove_artifacts(bp, trough_idxs, t)

# #peak_vals, new_peak_idxs = remove_random(bp, peak_idxs, t)

# manual_peak_detection(bp, t, peak_idxs, peak_vals)


# # fig, ax = plt.subplots()
# # ax.plot(t, bp, color="orange")
# # points = ax.scatter(t[peak_idxs], bp[peak_idxs], color="blue")
# # #plt.scatter(t[trough_idxs], trough_vals, color="green", s=10)
# # ax.set_title("BP with plotted Systolic Peaks random outliers filtered out")
# # plt.show()

# print("times after", t[peak_idxs])
# print("peaks after", bp[peak_idxs])
# peak_vals, peak_times = interpolate_missed_beats(bp, t, bp[peak_idxs], peak_idxs, 20, 1.2)

# # peak_vals2, peak_times2 = interpolate_missed_beats_part2(new_peak_idxs, bp[new_peak_idxs], bp, t)

# # print("test", peak_times1)
# # print("test", peak_times2)

# # #print(peak_vals)
# # #print(peak_times)


# peak_vals_uniform, t_uniform = interpolate_and_resample(peak_vals, peak_times, 2)
# # peak_vals_uniform2, t_uniform2 = interpolate_and_resample(peak_times2, peak_vals2)

# # print("peak vals uniform 1", peak_vals_uniform1)
# # print("peak vals uniform 2", peak_vals_uniform2)

# # plt.figure(figsize=(12, 6))
# # # Original signal
# # plt.subplot(3, 1, 1)
# # plt.plot(t, bp, label='Signal')
# # plt.scatter(t_uniform1, peak_vals_uniform1)
# # plt.title('one beat')
# # plt.xlabel('Time (s)')
# # plt.ylabel('Amplitude')
# # plt.grid()
# # plt.legend()

# #     # Filtered signal
# # plt.subplot(3, 1, 2)
# # plt.plot(t, bp, label='Signal', color='orange')
# # plt.scatter(t_uniform2, peak_vals_uniform2)
# # plt.title('multiple interpolate beats')
# # plt.xlabel('Time (s)')
# # plt.ylabel('Amplitude')
# # plt.grid()
# # plt.legend()

# # plt.tight_layout()
# # plt.show()

# #trough_vals_uniform, t_uniform_t = interpolate_and_resample(t[trough_idxs], trough_vals)

# plt.plot(t, bp, color="orange")
# plt.scatter(t_uniform, peak_vals_uniform, color="blue")
# #plt.scatter(t_uniform_t, trough_vals_uniform, color="green", s=10)
# plt.title("BP with plotted Systolics outliers filtered out + interpolation")
# plt.show()

