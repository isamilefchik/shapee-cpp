/**
 * @author      : Isa Milefchik (isavmilefchik@gmail.com)
 * @created     : Saturday Jun 06, 2020 23:19:36 PDT
 * @file        : shapee
 */

#include "shapee.hpp"
#include <iostream>
#include <mutex>
#include "/opt/homebrew/Cellar/libomp/11.0.1/include/omp.h"

Shapee::Shapee(int window_size, int hop_size, float w)
    : _window_size(window_size), _hop_size(hop_size), _w(w)
{
    _fftw_double_data = (double*) fftw_malloc(sizeof(double) * _window_size);
    _fftw_complex_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * _window_size);
    _fft_plan = fftw_plan_dft_r2c_1d(_window_size, _fftw_double_data, _fftw_complex_data, FFTW_ESTIMATE);
    _ifft_plan = fftw_plan_dft_c2r_1d(_window_size, _fftw_complex_data, _fftw_double_data, FFTW_ESTIMATE);
}

Shapee::~Shapee()
{
    fftw_destroy_plan(_fft_plan);
    fftw_destroy_plan(_ifft_plan);
    fftw_free(_fftw_double_data);
    fftw_free(_fftw_complex_data);
}

Shapee::AudioBuffer Shapee::shape(AudioBuffer& freq_src, AudioBuffer& ampl_src)
{

    AudioBuffer result(freq_src.size());

    // Check that they have the same number of channels
    assert(freq_src.size() == ampl_src.size());

    // For each channel
    for (int c = 0; c < freq_src.size(); ++c) {

        // Calculate freq STFT
        std::cout << "\tChannel " << c << ": Calculating frequency STFT...\t"
                  << std::flush;
        Spectrogram freq_spec = compute_stft(freq_src.at(c));
        std::cout << "done." << std::endl;

        // Calculate amplitude STFT
        std::cout << "\tChannel " << c << ": Calculating amplitude STFT...\t"
                  << std::flush;
        Spectrogram ampl_spec = compute_stft(ampl_src.at(c));
        std::cout << "done." << std::endl;

        // Begin frequency shaping process
        std::cout << "\tChannel " << c << ": Frequency shaping...\t\t\t"
                  << std::flush;
        int N = freq_spec.at(0).size();

        int result_len = std::min(freq_spec.size(), ampl_spec.size());
        Spectrogram result_spec(result_len);

        std::mutex result_lock;

#pragma omp parallel for
        for (int i = 0; i < result_len; ++i) {
            // Compute polar form of spectrogram frame
            FFTBins freq_polar(freq_spec.at(i).size());
            FFTBins ampl_polar(ampl_spec.at(i).size());
            for (size_t j = 0; j < freq_spec.at(i).size(); ++j) {
                freq_polar.at(j) = complex_to_polar(freq_spec.at(i).at(j));
                ampl_polar.at(j) = complex_to_polar(ampl_spec.at(i).at(j));
            }

            // Equation 4.6: Compute S(j)
            // ==========================

            int S_range = ((N / 2) - 1) / _w;
            std::vector<float> S(S_range + 1);

            for (int j = 0; j <= S_range; ++j) {
                float f_sum = 0, a_sum = 0;

                for (int n = 0; n <= _w; ++n) {
                    f_sum += freq_polar.at(static_cast<int>(j * _w) + n).real();
                    a_sum += ampl_polar.at(static_cast<int>(j * _w) + n).real();
                }

                S.at(j) = a_sum / f_sum;
            }

            // Equation 4.7 & 4.8: Compute freq-shaped STFT
            // ============================================
            FFTBins result_polar(N);

            for (int k = 0; k <= N / 2 - 1; ++k) {
                int k_prime = k / _w;

                float r_real = freq_polar.at(k).real() * S.at(k_prime);
                float r_imag = freq_polar.at(k).imag();

                result_polar.at(k) = std::complex<float>(r_real, r_imag);

                // Reflect on upper-half
                result_polar.at(N - 1 - k) =
                    std::complex<float>(r_real, r_imag);
            }

            // Convert and copy complex form into result spectrogram
            result_spec.at(i).resize(N);
            for (size_t j = 0; j < result_polar.size(); ++j) {
                result_spec.at(i).at(j) = polar_to_complex(result_polar.at(j));
            }
        }

        std::cout << "done." << std::endl;

        std::cout << "\tChannel " << c << ": Calculating shaped iSTFT...\t\t"
                  << std::flush;
        Wave result_wave = compute_istft(result_spec);
        std::cout << "done." << std::endl;

        result_lock.lock();
        result.at(c) = median_filter(result_wave, 3);
        result_lock.unlock();
    }

    return normalize_audio(result);
}

Shapee::Spectrogram Shapee::compute_stft(Wave& wave)
{
    // Compute number of windows needed
    int num_windows = std::ceil(wave.size() / static_cast<float>(_hop_size));

    // Reserve memory for spectrogram
    Spectrogram spec;
    spec.resize(num_windows);

    std::mutex spec_lock;

// #pragma omp parallel for
    for (int i = 0; i < wave.size(); i += _hop_size) {
        // Get window
        Wave window = get_window(wave, i);

        // Compute FFT of window
        FFTBins temp = compute_fft(window);

        spec_lock.lock();
        spec.at(i / _hop_size) = temp;
        spec_lock.unlock();
    }

    return spec;
}

Shapee::Wave Shapee::compute_istft(Spectrogram& spec)
{
    // Reserve memory for wave signal
    int wave_size = (spec.size() * _hop_size) + (_window_size - _hop_size);
    Wave wave(wave_size, 0.f);

    std::mutex wave_lock;

// #pragma omp parallel for
    for (int i = 0; i < spec.size(); ++i) {
        // Compute inverse FFT of current frame
        Wave window = compute_ifft(spec.at(i));

        // Add to output wave
        wave_lock.lock();
        int base_idx = _hop_size * i;
        for (int j = 0; j < window.size(); ++j) {
            wave.at(base_idx + j) += window.at(j);
        }
        wave_lock.unlock();
    }

    return wave;
}

// Hamming window function
Shapee::Wave Shapee::get_window(const Wave& wave, int base_idx)
{
    Wave window(_window_size, 0);

    static const float scale_factor = 2 * M_PI / _window_size;
    static const float a_0 = 0.54f;
    static const float a_1 = 1.f - a_0;

    // Apply Hamming window
    for (int j = 0; j < _window_size; ++j) {
        float wave_sample =
            base_idx + j < wave.size() ? wave.at(base_idx + j) : 0;
        window.at(j) = (a_0 - (a_1 * std::cos(scale_factor * j))) * wave_sample;
    }

    return window;
}

Shapee::Wave Shapee::median_filter(const Wave& wave, int width) {
    Wave result(wave.size(), 0);

#pragma omp parallel for
    for (int i = 0; i < wave.size(); ++i) {
        int beg_idx = std::max(0, i - width);
        int end_idx = std::min(static_cast<int>(wave.size()-1), i + width);

        int n = end_idx - beg_idx + 1;
        float window[n];
        for (int j = 0; j < n; ++j) {
            window[j] = wave.at(beg_idx + j);
        }

        result.at(i) = calc_median(window, n);
    }

    return result;
}

float Shapee::calc_median(float window[], int n) {
    std::sort(window, window+n);

    // check for even case 
    if (n % 2 != 0) 
       return (float) window[n/2]; 
      
    return (float) (window[(n-1)/2] + window[n/2]) / 2.0; 
}

Shapee::AudioBuffer Shapee::normalize_audio(const AudioBuffer& audio) {

    // Calculate absolute max value:
    float max_val = 0;
    for (int c = 0; c < audio.size(); ++c) {
        for (int i = 0; i < audio.at(c).size(); ++i) {
            if (std::abs(audio.at(c).at(i)) > max_val) {
                max_val = std::abs(audio.at(c).at(i));
            }
        }
    }

    AudioBuffer normalized_audio(audio.size());

    // Normalize:
    for (int c = 0; c < audio.size(); ++c) {
        normalized_audio.at(c).resize(audio.at(c).size());
#pragma omp parallel for
        for (int i = 0; i < audio.at(c).size(); ++i) {
            normalized_audio.at(c).at(i) = std::clamp(audio.at(c).at(i) / max_val, -1.f, 1.f);
        }
    }

    return normalized_audio;
}

Shapee::FFTBins Shapee::compute_fft(Wave& wave)
{
    FFTBins bins(wave.size());

    for (int i = 0; i < wave.size(); ++i) {
        _fftw_double_data[i] = wave.at(i);
    }

    fftw_execute(_fft_plan);

    for (int i = 0; i < wave.size(); ++i) {
        bins.at(i) = std::complex<float>(_fftw_complex_data[i][0], _fftw_complex_data[i][1]);
    }

    return bins;
}

Shapee::Wave Shapee::compute_ifft(FFTBins& bins)
{
    Wave wave(bins.size());

    for (int i = 0; i < bins.size(); ++i) {
        _fftw_complex_data[i][0] = bins.at(i).real();
        _fftw_complex_data[i][1] = bins.at(i).imag();
    }

    fftw_execute(_ifft_plan);

    for (int i = 0; i < wave.size(); ++i) {
        wave.at(i) = _fftw_double_data[i];
    }

    return wave;
}

std::complex<float> Shapee::complex_to_polar(std::complex<float> c)
{
    return std::complex<float>(std::abs(c), std::arg(c));
}

std::complex<float> Shapee::polar_to_complex(std::complex<float> p)
{
    return std::polar(p.real(), p.imag());
}
