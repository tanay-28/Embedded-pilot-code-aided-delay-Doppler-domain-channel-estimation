# Embedded Pilot Code-Aided Delay-Doppler Domain Channel Estimation

MATLAB implementation of the embedded pilot code-aided channel estimation scheme for Orthogonal Time-Frequency Space (OTFS) modulation, as proposed in:

> **T. Agrawal and R. Kadlimatti**, "Embedded Pilot Code-Aided Delay-Doppler Domain Channel Estimation," *2026 IEEE 23rd Consumer Communications & Networking Conference (CCNC)*, Las Vegas, NV, 2026. DOI: 10.1109/CCNC65079.2026.11366337

---

## Overview

OTFS modulation transforms time-varying wireless channels into a sparse, flat-fading 2D delay-Doppler (DD) grid, making it well-suited for high-mobility environments where OFDM performs poorly. Accurate channel state information (CSI) is essential for reliable OTFS performance.

Conventional embedded pilot-aided channel estimation places a single amplified pilot symbol in the DD transmit grid. While spectrally efficient, this approach requires significantly higher pilot power to distinguish channel taps from the noise floor, resulting in an increased peak-to-average power ratio (PAPR) and susceptibility to channel fades at the pilot location.

This work proposes replacing the single pilot symbol with a **Barker code sequence** embedded in the DD transmit grid. At the receiver, cross-correlation with the known pilot code is used to estimate the DD channel transfer function. The ideal autocorrelation properties of the Barker code allow accurate tap detection without increasing pilot power, achieving:

- Lower PAPR (~6 dB reduction) compared to the single pilot scheme
- Accurate channel estimation at low SNRs
- Convergence to ideal channel estimation and BER performance

---

## System Parameters

| Parameter | Value |
|-----------|-------|
| Carrier frequency | 2 GHz |
| Subcarrier spacing | 15 kHz |
| Symbols per frame (N) | 16 |
| Subcarriers (M) | 128 |
| Bandwidth | 1.92 MHz |
| Modulation | 4-QAM |
| Channel taps | 4 |
| Max delay tap (lτ) | 3 |
| Max Doppler tap (kv) | 2 |
| Pilot code | Barker, length L = 11 |
| Channel model | Synthetic 3GPP |

---

## Method

### Transmitter
A Barker code sequence `xb` of length `L` is placed in the DD transmit grid at the pilot location `(kp, lp)` where `kp = N/2` and `lp = M/2`:

```
x[k,l] = xb,        k = kp, lp - lg ≤ l ≤ lp + lg
        = 0,         guard region around pilot
        = xd[k,l],   data symbols elsewhere
```

where `lg = (L-1)/2` is the guard length determined by the code length.

### Receiver
The cross-correlation of the received DD symbols with the known Barker code is computed:

```
R_y,x[ki, li + C] = sum_{n=-lg}^{lg} Y_DD[ki, li+n] * xb*[n]
```

Channel taps are detected by thresholding the cross-correlation output. When `|R_y,x| ≥ τ`, a path is declared present and the channel coefficient is estimated as `R_y,x / L`. The estimated channel is then used with a time-domain LMMSE detector for data equalization.

---

## Repository Structure

```
├── p3_coded_estimation_main.m          # Main simulation script (P3 coded pilot)
├── barker_coded_pilot_estimation.m     # Cross-correlation based channel estimation
├── embedded_pilot_channel_estimation.m # Single pilot (baseline) channel estimation
├── Synthetic_channel_gen.m             # 3GPP synthetic channel generation
├── Block_LMMSE_detector.m              # Time-domain LMMSE equalizer
├── otfs_modulation.m                   # OTFS modulator (DD → time domain)
├── otfs_demodulation.m                 # OTFS demodulator (time domain → DD)
├── delay_time_transmission.m           # Delay-time domain channel application
├── Gen_time_domain_channel.m           # Time-domain channel matrix generation
├── Generate_2D_data_grid.m             # DD grid data placement
├── transmit_data_gen.m                 # QAM data symbol generation
├── papr_calculation.m                  # PAPR computation
└── p3.m                                # P3 polyphase code generation
```

---

## Running the Simulation

1. Clone the repository and open MATLAB
2. Navigate to this folder
3. Run the main script:

```matlab
p3_coded_estimation_main
```

The script will generate three figures:
- **BER vs SNR** — comparing uncoded pilot, Barker coded pilot, and ideal channel estimation
- **PAPR CCDF** — comparing peak power distribution between schemes
- **Cross-correlation output** — visualizing the DD domain channel tap detection

---

## Results

The Barker code-based pilot scheme achieves:

- **BER convergence to ideal** channel estimation at matched total transmit power
- **~6 dB PAPR reduction** compared to the single amplified pilot symbol scheme
- **Lower RMSE** in channel coefficient estimation across all SNR levels

---

## Citation

If you use this code, please cite:

```bibtex
@inproceedings{agrawal2026embedded,
  title={Embedded Pilot Code-Aided Delay-Doppler Domain Channel Estimation},
  author={Agrawal, Tanay and Kadlimatti, Ravi},
  booktitle={2026 IEEE 23rd Consumer Communications \& Networking Conference (CCNC)},
  year={2026},
  address={Las Vegas, NV},
  doi={10.1109/CCNC65079.2026.11366337}
}
```

---

## References

[1] P. Raviteja, K. T. Phan and Y. Hong, "Embedded Pilot-Aided Channel Estimation for OTFS in Delay-Doppler Channels," *IEEE Transactions on Vehicular Technology*, vol. 68, no. 5, pp. 4906-4917, 2019.

[2] R. Hadani et al., "Orthogonal Time Frequency Space Modulation," *2017 IEEE WCNC*, San Francisco, CA, 2017.

[3] N. Levanon and E. Mozeson, "Phase-Coded Pulse," in *Radar Signals*, Chapter 6, John Wiley & Sons, 2004.
