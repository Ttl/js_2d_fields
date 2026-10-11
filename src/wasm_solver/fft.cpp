// Real FFTs for the time-domain (eye) simulation.
//
// Eigen's unsupported FFT module with its default header-only kissfft backend.
// Transforms of one size share the plan Eigen::FFT caches for it.
//
// Spectra are half spectra of n/2 + 1 interleaved (re, im) doubles. n must be even.
// The inverse includes the 1/n factor (x[k] = 1/n sum X[m] e^{+2 pi i k m / n}).

#include <unsupported/Eigen/FFT>
#include <complex>
#include <emscripten.h>

static Eigen::FFT<double>& engine() {
    static Eigen::FFT<double> fft;
    static bool init = false;
    if (!init) { fft.SetFlag(Eigen::FFT<double>::HalfSpectrum); init = true; }
    return fft;
}

extern "C" {

// nbatch inverse transforms of size n: spec holds nbatch half spectra back to back,
// out receives nbatch real signals of n samples. Returns 0, or -1 for a bad size.
EMSCRIPTEN_KEEPALIVE
int irfft_batch(int n, int nbatch, const double* spec, double* out) {
    if (n < 2 || (n & 1) || nbatch < 1) return -1;
    auto& fft = engine();
    const int nh = n / 2 + 1;
    for (int b = 0; b < nbatch; b++) {
        const auto* src = reinterpret_cast<const std::complex<double>*>(spec + 2 * (size_t)b * nh);
        fft.inv(out + (size_t)b * n, src, n);
    }
    return 0;
}

// Forward transform of n real samples into n/2 + 1 complex bins (no scaling).
EMSCRIPTEN_KEEPALIVE
int rfft(int n, const double* in, double* spec) {
    if (n < 2 || (n & 1)) return -1;
    auto& fft = engine();
    fft.fwd(reinterpret_cast<std::complex<double>*>(spec), in, n);
    return 0;
}

}
