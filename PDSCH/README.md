# 5G NR PDSCH Simulation — EE5770 Course Project

## Overview
Complete 5G NR PDSCH link-level simulator comparing LDPC and convolutional codes
over L-tap frequency-selective Rayleigh fading with ZF/MMSE equalization.

## Directory Structure
```
PDSCH/
├── pdsch_sim.m              # Main simulation (run this)
├── report.tex               # LaTeX report
├── README.md                # This file
├── encoder/                 # Encoding modules
│   ├── crc_utils.m          # CRC24A/CRC24B encode/check
│   ├── code_block_segment.m # Code block segmentation
│   ├── ldpc_encode_block.m  # 5G NR LDPC encoder (BG1/BG2)
│   └── conv_encode.m        # Rate-1/2 convolutional encoder
├── decoder/                 # Decoding modules
│   ├── ldpc_decode_block.m  # Min-Sum LDPC decoder
│   └── viterbi_decode.m     # Soft-decision Viterbi decoder
├── ofdm/                    # OFDM processing
│   ├── ofdm_modulate.m      # IFFT + CP insertion
│   └── ofdm_demod_equalize.m # CP removal + FFT + equalization
├── channel/                 # Channel model
│   └── channel_multipath.m  # L-tap Rayleigh fading + AWGN
├── equalizer/               # (Integrated in ofdm_demod_equalize.m)
│   └── README.md            # ZF and MMSE share common interface
├── simulation/              # Support modules
│   ├── rate_matching.m      # Rate match / rate recover
│   ├── scrambling.m         # Gold sequence scrambling
│   └── modulation.m         # Gray mapping / LLR demapping
└── tests/                   # Unit tests
    └── test_all_blocks.m    # All processing block tests
```

## Quick Start

### Run Full Simulation
```matlab
>> pdsch_sim
```
This runs all 5 experiments and generates all plots automatically:
- **Fig 1**: BER/BLER for Conv K=7, K=9, LDPC × ZF/MMSE (L=4)
- **Fig 2**: ZF vs MMSE at L=4 and L=16
- **Fig 3**: Per-subcarrier SNR CDF at Eb/N0 = 5 dB
- **Fig 4**: LLR mismatch (per-SC vs global SNR)
- **Fig 5**: Decode time comparison

### Run Unit Tests
```matlab
>> cd tests
>> test_all_blocks
```

### Save Plots for Report
After running `pdsch_sim`, save each figure:
```matlab
>> saveas(figure(1), 'fig1_ber_bler.png')
>> saveas(figure(2), 'fig2_zf_mmse.png')
>> saveas(figure(3), 'fig3_snr_cdf.png')
>> saveas(figure(4), 'fig4_llr_mismatch.png')
>> saveas(figure(5), 'fig5_complexity.png')
```

### Compile Report
```bash
pdflatex report.tex
```

## Configuration Parameters
All parameters are set in the `cfg` struct at the top of `pdsch_sim.m`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cfg.A` | 1000 | Transport block size (bits) |
| `cfg.R` | 1/2 | Target code rate |
| `cfg.M` | 2 | Bits/symbol: 2=QPSK, 4=16QAM, 6=64QAM |
| `cfg.NFFT` | 256 | FFT size |
| `cfg.L` | 4 | Number of channel taps |
| `cfg.EbN0_dB` | 0:1:8 | Eb/N0 sweep range (dB) |
| `cfg.num_frames` | 500 | Monte Carlo frames per SNR point |
| `cfg.ldpc_max_iter` | 8 | Max LDPC BP iterations |
| `cfg.conv_K` | 7 | Constraint length (7 or 9) |
| `cfg.conv_gen` | [133 171] | Generator polynomials (octal) |
| `cfg.coding` | 'conv' | Coding scheme: 'conv' or 'ldpc' |
| `cfg.equalizer` | 'mmse' | Equalizer: 'mmse' or 'zf' |
| `cfg.use_global_snr` | false | LLR mismatch experiment flag |

## Coding Schemes
- **Scheme A — 5G NR LDPC**: BG1/BG2, lifted parity-check matrix, Min-Sum decoder (α=0.75)
- **Scheme B — Conv K=7**: (133,171)₈, 64-state trellis, soft Viterbi
- **Scheme B — Conv K=9**: (561,753)₈, 256-state trellis, soft Viterbi

## Equalizer Interface
Both equalizers accept `(Y_k, H_k, σ²_w)` and return `(Ŝ_k, γ_k)`:
- **ZF**: `Ŝ = Y/H`, `γ = |H|²/σ²`
- **MMSE**: `Ŝ = H*Y/(|H|²+σ²)`, `γ = |H|⁴/(σ²(|H|²+σ²))`

## Requirements
- MATLAB R2020a or later
- No additional toolboxes required
