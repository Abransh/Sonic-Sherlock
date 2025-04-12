#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "../fingerprinter/fft.h"
#include "../fingerprinter/peaks.h"
#include "../fingerprinter/fingerprint.h"

namespace py = pybind11;

// Helper function to convert numpy array to std::vector
template<typename T>
std::vector<T> numpy_to_vector(py::array_t<T> array) {
    auto buffer = array.request();
    T* ptr = static_cast<T*>(buffer.ptr);
    return std::vector<T>(ptr, ptr + buffer.size);
}

// Helper to convert 2D numpy array to std::vector<std::vector<double>>
std::vector<std::vector<double>> numpy_to_2d_vector(py::array_t<double> array) {
    auto buffer = array.request();
    if (buffer.ndim != 2) {
        throw std::runtime_error("Number of dimensions must be 2");
    }
    
    std::vector<std::vector<double>> result(buffer.shape[0], std::vector<double>(buffer.shape[1]));
    double* ptr = static_cast<double*>(buffer.ptr);
    
    for (size_t i = 0; i < buffer.shape[0]; ++i) {
        for (size_t j = 0; j < buffer.shape[1]; ++j) {
            result[i][j] = ptr[i * buffer.shape[1] + j];
        }
    }
    
    return result;
}

// Helper to convert std::vector<std::vector<double>> to 2D numpy array
py::array_t<double> vector_to_numpy_2d(const std::vector<std::vector<double>>& data) {
    if (data.empty()) {
        return py::array_t<double>(0);
    }
    
    size_t rows = data.size();
    size_t cols = data[0].size();
    
    py::array_t<double> result({rows, cols});
    auto buffer = result.request();
    double* ptr = static_cast<double*>(buffer.ptr);
    
    for (size_t i = 0; i < rows; ++i) {
        if (data[i].size() != cols) {
            throw std::runtime_error("Input vector is not rectangular");
        }
        for (size_t j = 0; j < cols; ++j) {
            ptr[i * cols + j] = data[i][j];
        }
    }
    
    return result;
}

PYBIND11_MODULE(audiofp_core, m) {
    m.doc() = "Audio fingerprinting core C++ implementation";
    
    // Expose FFT class
    py::class_<audiofp::FFT>(m, "FFT")
        .def_static("compute", [](py::array_t<double> signal) {
            auto vec_signal = numpy_to_vector<double>(signal);
            auto result = audiofp::FFT::compute(vec_signal);
            
            // Convert complex vector to numpy array
            size_t n = result.size();
            py::array_t<std::complex<double>> numpy_result(n);
            auto buffer = numpy_result.request();
            std::complex<double>* ptr = static_cast<std::complex<double>*>(buffer.ptr);
            
            for (size_t i = 0; i < n; ++i) {
                ptr[i] = result[i];
            }
            
            return numpy_result;
        })
        .def_static("magnitude_spectrum", [](py::array_t<std::complex<double>> fft_result) {
            auto buffer = fft_result.request();
            std::complex<double>* ptr = static_cast<std::complex<double>*>(buffer.ptr);
            std::vector<std::complex<double>> vec_fft(ptr, ptr + buffer.size);
            
            auto result = audiofp::FFT::magnitude_spectrum(vec_fft);
            
            // Convert to numpy array
            py::array_t<double> numpy_result(result.size());
            auto result_buffer = numpy_result.request();
            double* result_ptr = static_cast<double*>(result_buffer.ptr);
            
            for (size_t i = 0; i < result.size(); ++i) {
                result_ptr[i] = result[i];
            }
            
            return numpy_result;
        })
        .def_static("apply_window", [](py::array_t<double> signal, int window_type) {
            auto vec_signal = numpy_to_vector<double>(signal);
            auto result = audiofp::FFT::apply_window(vec_signal, window_type);
            
            // Convert to numpy array
            py::array_t<double> numpy_result(result.size());
            auto buffer = numpy_result.request();
            double* ptr = static_cast<double*>(buffer.ptr);
            
            for (size_t i = 0; i < result.size(); ++i) {
                ptr[i] = result[i];
            }
            
            return numpy_result;
        });
    
    // Expose Peak struct
    py::class_<audiofp::Peak>(m, "Peak")
        .def(py::init<int, int, double>())
        .def_readwrite("time", &audiofp::Peak::time)
        .def_readwrite("frequency", &audiofp::Peak::frequency)
        .def_readwrite("magnitude", &audiofp::Peak::magnitude);
    
    // Expose PeakFinder class
    py::class_<audiofp::PeakFinder>(m, "PeakFinder")
        .def_static("find_peaks", [](py::array_t<double> spectrogram, int neighborhood_size, double threshold, int max_peaks) {
            auto vec_spectrogram = numpy_to_2d_vector(spectrogram);
            auto result = audiofp::PeakFinder::find_peaks(vec_spectrogram, neighborhood_size, threshold, max_peaks);
            return result;
        })
        .def_static("generate_constellation", &audiofp::PeakFinder::generate_constellation);
    
    // Expose FingerprintHash struct
    py::class_<audiofp::FingerprintHash>(m, "FingerprintHash")
        .def(py::init<uint32_t, uint32_t>())
        .def_readwrite("hash", &audiofp::FingerprintHash::hash)
        .def_readwrite("time_offset", &audiofp::FingerprintHash::time_offset);
    
    // Expose Fingerprinter class
    py::class_<audiofp::Fingerprinter>(m, "Fingerprinter")
        .def_static("generate_spectrogram", [](py::array_t<double> samples, int sample_rate, int window_size, int hop_size) {
            auto vec_samples = numpy_to_vector<double>(samples);
            auto result = audiofp::Fingerprinter::generate_spectrogram(vec_samples, sample_rate, window_size, hop_size);
            return vector_to_numpy_2d(result);
        })
        .def_static("generate_fingerprints", [](py::array_t<double> spectrogram, int sample_rate) {
            auto vec_spectrogram = numpy_to_2d_vector(spectrogram);
            auto result = audiofp::Fingerprinter::generate_fingerprints(vec_spectrogram, sample_rate);
            return result;
        })
        .def_static("fingerprint_audio", [](py::array_t<double> samples, int sample_rate) {
            auto vec_samples = numpy_to_vector<double>(samples);
            auto result = audiofp::Fingerprinter::fingerprint_audio(vec_samples, sample_rate);
            return result;
        })
        .def_static("match_fingerprints", &audiofp::Fingerprinter::match_fingerprints);
}