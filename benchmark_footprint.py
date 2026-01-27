import numpy as np
import scipy.stats as stats
import time

def original_footprint_correction(alpha_i, beam_size, sample_size):
    samplesize_microns = sample_size * 10000

    start_time = time.time()
    footprint = np.array([
        0.5 * beam_size / np.sin(np.deg2rad(alpha)) if alpha != 0 else 0.5 * beam_size / np.sin(np.deg2rad(1e-3))
        for alpha in alpha_i
    ])

    beam_fraction = np.array([
        (stats.norm.cdf(samplesize_microns / 2, 0, ftp) - stats.norm.cdf(-samplesize_microns / 2, 0, ftp))
        for ftp in footprint
    ])
    end_time = time.time()
    return end_time - start_time, beam_fraction

def optimized_footprint_correction(alpha_i, beam_size, sample_size):
    samplesize_microns = sample_size * 10000

    start_time = time.time()

    # Vectorized footprint calculation
    alpha_rad = np.deg2rad(alpha_i)

    # Handle alpha_i == 0
    # Use np.where to select the divisor angle
    # The original code logic: if alpha != 0 -> alpha, else -> 1e-3
    effective_alpha_rad = np.where(alpha_i != 0, alpha_rad, np.deg2rad(1e-3))

    footprint = 0.5 * beam_size / np.sin(effective_alpha_rad)

    # Vectorized beam_fraction
    # stats.norm.cdf works with arrays. scale=footprint
    beam_fraction = (stats.norm.cdf(samplesize_microns / 2, 0, footprint) - stats.norm.cdf(-samplesize_microns / 2, 0, footprint))

    end_time = time.time()
    return end_time - start_time, beam_fraction

# Setup data
# Using a reasonably large size to see performance differences
N = 10000
alpha_i = np.linspace(0, 10, N) # 0 to 10 degrees, includes 0
beam_size = 9.6
sample_size = 1.0

print(f"Benchmarking with N={N} elements...")

# Warmup / Run Original
t_orig, res_orig = original_footprint_correction(alpha_i, beam_size, sample_size)
print(f"Original Time: {t_orig:.6f}s")

# Warmup / Run Optimized
t_opt, res_opt = optimized_footprint_correction(alpha_i, beam_size, sample_size)
print(f"Optimized Time: {t_opt:.6f}s")

# Validation
is_close = np.allclose(res_orig, res_opt)
print(f"Match: {is_close}")
if not is_close:
    print("Max difference:", np.max(np.abs(res_orig - res_opt)))

if t_opt > 0:
    print(f"Speedup: {t_orig / t_opt:.2f}x")
else:
    print("Optimized time is 0")
