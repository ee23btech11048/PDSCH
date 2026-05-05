% Equalizer Module
% The ZF and MMSE equalizers are integrated in ofdm/ofdm_demod_equalize.m
% with a common interface:
%
%   [eq_sym, snr_out] = ofdm_demod_equalize(rx_sig, h_taps, NFFT, NCP, 
%                                            num_sym, sc_idx, sigma2_w, eq_type)
%
% Inputs:
%   rx_sig   - received time-domain signal
%   h_taps   - channel impulse response (L x 1)
%   NFFT     - FFT size
%   NCP      - cyclic prefix length
%   num_sym  - number of OFDM symbols
%   sc_idx   - data subcarrier indices
%   sigma2_w - noise variance
%   eq_type  - 'zf' or 'mmse'
%
% Outputs:
%   eq_sym   - equalized symbols
%   snr_out  - per-subcarrier post-equalization SNR
