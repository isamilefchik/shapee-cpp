/**
 * @author      : Isa Milefchik (isavmilefchik@gmail.com)
 * @created     : Saturday Jun 06, 2020 23:14:46 PDT
 * @file        : main
 */

#include <iostream>
#include "audiofile.hpp"
#include "shapee.hpp"
#include "cxxopts.hpp"

int main(int argc, char** argv)
{
    /*
     * ===========================================================================
     * ARGUMENT PARSING
     * ===========================================================================
     */

    std::cout << "\nWelcome to...\n\n";
    std::string banner =
    " |===================================================| \n"
    " |  _____ _____ _____ _____ _____ _____   _     _    | \n"
    " | |   __|  |  |  _  |  _  |   __|   __|_| |_ _| |_  | \n"
    " | |__   |     |     |   __|   __|   __|_   _|_   _| | \n"
    " | |_____|__|__|__|__|__|  |_____|_____| |_|   |_|   | \n"
    " |                                                   | \n"
    " |===================================================| \n";
    std::cout << banner << std::endl;


    cxxopts::Options options("shapee++", "An audio shaping and synthesis tool.");

    options.add_options()
        (
         "f,freq_src",
         "Path to audio file for frequency component of synthesis.",
         cxxopts::value<std::string>()
        )
        (
         "t,timbre_src",
         "Path to audio file for timbre component of synthesis.",
         cxxopts::value<std::string>()
        )
        (
         "o,output",
         "Path to output file.",
         cxxopts::value<std::string>()->default_value("./shapee_output.wav")
        )
        (
         "w,omega",
         "Value for omega in frequency shaping formulation.",
         cxxopts::value<float>()->default_value("5")
        )
        (
         "stft_len",
         "Length of STFT windows (e.g., 256, 512, 1024, etc.).",
         cxxopts::value<int>()->default_value("2048")
        )
        (
         "stft_hop",
         "Hop size of STFT windows (e.g., 256, 512, 1024, etc.).",
         cxxopts::value<int>()->default_value("256")
        )
        ("h,help", "Print help message.")
        ;

    auto options_result = options.parse(argc, argv);
    if (options_result.count("help")) {
      std::cout << options.help() << std::endl;
      exit(0);
    }

    if (!options_result.count("freq_src")) {
        std::cout << "Please provide a path to the frequency audio source.\n\n";
        std::cout << options.help() << std::endl;
        exit(0);
    }

    if (!options_result.count("timbre_src")) {
        std::cout << "Please provide a path to the timbre audio source.\n\n";
        std::cout << options.help() << std::endl;
        exit(0);
    }

    std::string freq_src_path = options_result["freq_src"].as<std::string>();
    std::string timbre_src_path = options_result["timbre_src"].as<std::string>();
    std::string output_path = options_result["output"].as<std::string>();
    float omega = options_result["omega"].as<float>();
    int stft_len = options_result["stft_len"].as<int>();
    int stft_hop = options_result["stft_hop"].as<int>();

    /*
     * ===========================================================================
     * PROGRAM
     * ===========================================================================
     */

    AudioFile<float> af_freq, af_timbre;

    // Load freq audio
    std::cout << "Loading " << freq_src_path << " into memory...\t" << std::flush;
    af_freq.load(freq_src_path);
    std::cout << "done.\n";

    // Load timbre audio
    std::cout << "Loading " << timbre_src_path << " into memory...\t" << std::flush;
    af_timbre.load(timbre_src_path);
    std::cout << "done.\n\n";

    // Get audio format info
    int source_bit_depth = af_freq.getBitDepth();
    int source_num_channels = af_freq.getNumChannels();
    int source_sample_rate = af_freq.getSampleRate();

    // Frequency shape
    std::cout << "\U0001f608 Shaping...\n\n";
    Shapee my_shapee(stft_len, stft_hop, omega);
    Shapee::AudioBuffer buff = my_shapee.shape(af_freq.samples, af_timbre.samples);
    std::cout << "\nDone shaping.\n\n";

    // Write result
    std::cout << "Writing to " << output_path << "... " << std::flush;
    AudioFile<float> result;
    result.setBitDepth(source_bit_depth);
    result.setNumChannels(source_num_channels);
    result.setSampleRate(source_sample_rate);
    result.setAudioBuffer(buff);
    result.save(output_path);
    std::cout << "done.\n\n";

    return 0;
}

