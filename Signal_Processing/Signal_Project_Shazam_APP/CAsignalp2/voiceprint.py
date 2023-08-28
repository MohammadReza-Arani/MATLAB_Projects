from pydub import AudioSegment as audio
import numpy as np
from scipy import signal
from matplotlib import pyplot as plt


def audioread(file_name):
    sound = audio.from_mp3(file_name)
    # raw_data = sound.raw_data
    sample_rate = sound.frame_rate
    channels = sound.channels
    if channels == 1:
        clip = np.array(sound.get_array_of_samples())
    elif channels == 2:
        sound = sound.split_to_mono()
        channel_one = np.array(sound[0].get_array_of_samples())
        channel_two = np.array(sound[0].get_array_of_samples())
        clip = np.array([channel_one, channel_two])
    return clip, sample_rate


def voiceprint(clip, fs):

    return peaks


def peak_to_pair(peaks):
    delta_t = 20    # Time interval after the peak. (in pixels)
    delta_f = 10    # Distance from the frequency fo the peak. (in pixels)
    plot_lines = True
    time_lower_band = 0
    fanout = 3
    time_upper_band = time_lower_band + delta_t

    peaks = peaks[:256, :]
    t, f = np.nonzero(np.transpose(peaks))
    num_of_peaks = f.size

    pairs = np.ndarray(shape=(0, 4))

    for p in range(num_of_peaks):
        current_peak_f = f[p]
        current_peak_t = t[p]
        peak_pair_position = t != t[p]
        peak_pair_position = np.logical_and(
            peak_pair_position, (t > current_peak_t + time_lower_band))
        peak_pair_position = np.logical_and(
            peak_pair_position, (t < current_peak_t + time_upper_band))
        peak_pair_position = np.logical_and(
            peak_pair_position, (f > current_peak_f - delta_f))
        peak_pair_position = np.logical_and(
            peak_pair_position, (f < current_peak_f + delta_f))
        peak_pair_position = np.nonzero(peak_pair_position.astype(np.int))

        if peak_pair_position[0].size > fanout:
            peak_pair_position = peak_pair_position[0][:fanout]
        else:
            peak_pair_position = peak_pair_position[0]
        for i in peak_pair_position:
            feature_tuple = np.array([t[p], t[i], f[p], f[i]])
            pairs = np.vstack((pairs, feature_tuple))

    if plot_lines:
        fig, ax = plt.subplots()
        ax.imshow(peaks, cmap=plt.cm.binary)
        for p in pairs:
            ax.plot(p[:2], p[2:4])
        ax.set_title('Spectrogram local peaks and peak pairs')
        ax.set_xlabel('Time')
        ax.set_ylabel('Frequency')
        plt.show()


if __name__ == '__main__':
    clip, fs = audioread('viva.mp3')
    peaks = voiceprint(clip, fs)
    peak_to_pair(peaks)
