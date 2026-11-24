import matlab.engine
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent
MATLAB_DIR = PROJECT_DIR / "matlab functions"



def get_lf_powers(signal, fs):
    fs = float(fs)
    signal_matlab_array = matlab.double(signal.tolist())
    eng = matlab.engine.start_matlab()
    eng.addpath(str(MATLAB_DIR), nargout=0)
    lf_powers, entire_lf_within_coi = eng.get_lf_powers(signal_matlab_array, fs, nargout=2)
    return lf_powers, entire_lf_within_coi

def get_avg_lf(lf_powers, sample_seconds, fs, entire_lf_within_coi, start_idx, end_idx):
    fs = float(fs)
    sample_seconds = float(sample_seconds)

    eng = matlab.engine.start_matlab()
    eng.addpath(str(MATLAB_DIR), nargout=0)
    avg_lf = eng.get_avg_lf(lf_powers, sample_seconds, fs, entire_lf_within_coi, start_idx, end_idx, nargout=1)
    return avg_lf  
