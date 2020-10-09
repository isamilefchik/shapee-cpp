/**
 * @author      : Isa Milefchik (isavmilefchik@gmail.com)
 * @created     : Saturday Jun 06, 2020 23:19:53 PDT
 * @file        : shapee
 */

#ifndef SHAPEE_HPP
#define SHAPEE_HPP

#include <complex>
#include <vector>

class Shapee
{
  public:
    typedef std::vector<float> Wave;
    typedef std::vector<Wave> AudioBuffer;
    typedef std::vector<std::complex<float>> FFTBins;
    typedef std::vector<FFTBins> Spectrogram;

    Shapee(int window_size, int hop_size, float w);

    AudioBuffer shape(AudioBuffer& freq_src, AudioBuffer& ampl_src);

  private:
    Spectrogram compute_stft(Wave& wave);
    Wave compute_istft(Spectrogram& spec);
    FFTBins compute_fft(Wave& wave);
    // FFTBins compute_complex_fft(FFTBins wave);
    Wave compute_ifft(FFTBins& bins);
    Wave get_window(const Wave& wave, int i);
    Wave median_filter(const Wave& wave, int width);
    float calc_median(float window[], int n);
    AudioBuffer normalize_audio(const AudioBuffer& audio);
    std::complex<float> complex_to_polar(std::complex<float> c);
    std::complex<float> polar_to_complex(std::complex<float> p);

    int _window_size, _hop_size;
    float _w;
};

#endif /* SHAPEE_HPP */
