# shapee++

My original shapee program was too slow. So I rewrote it in C++ and used the Intel Math Kernel Library (MKL) for FFT calculations.

### Usage:

```
Welcome to...

 |===================================================|
 |  _____ _____ _____ _____ _____ _____   _     _    |
 | |   __|  |  |  _  |  _  |   __|   __|_| |_ _| |_  |
 | |__   |     |     |   __|   __|   __|_   _|_   _| |
 | |_____|__|__|__|__|__|  |_____|_____| |_|   |_|   |
 |                                                   |
 |===================================================|

An audio shaping and synthesis tool.
Usage:
  shapee++ [OPTION...]

  -f, --freq_src arg    Path to audio file for frequency component of
                        synthesis.
  -t, --timbre_src arg  Path to audio file for timbre component of synthesis.
  -o, --output arg      Path to output file. (default: ./shapee_output.wav)
  -w, --omega arg       Value for omega in frequency shaping formulation.
                        (default: 5)
      --stft_len arg    Length of STFT windows (e.g., 256, 512, 1024, etc.).
                        (default: 2048)
      --stft_hop arg    Hop size of STFT windows (e.g., 256, 512, 1024,
                        etc.). (default: 256)
  -h, --help            Print help message.
```

### References:

1. C. Penrose, "Frequency Shaping of Audio Signals" (ICMC 2001) [link](https://quod.lib.umich.edu/cgi/p/pod/dod-idx/frequency-shaping-of-audio-signals.pdf?c=icmc&format=pdf&idno=bbp2372.2001.082)
