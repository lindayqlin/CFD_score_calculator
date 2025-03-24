# Cutting frequency determination (CFD) score calculator
- A simple, standalone CFD score calculator that supports gaps in RNA-DNA alignment
- Implemented in both Rust (faster) and Python (original)
- CFD score is an empirical measure of the likelihood of a given guide RNA (spacer) sequence complexed with SpCas9 cleaving a given DNA sequence. CFD scores are helpful for prioritizing candidate off-target sites for further evaluation. They range from 0 (minimal cutting potential) to 1 (for on-target site)

## Usage
### Rust
```rust
// Import:
mod cfd_score_calculator;
// Usage:
// cfd_score_calculator::calculate_cfd(spacer, protospacer, pam);
// Example:
cfd_score_calculator::calculate_cfd("CTAACAGTTGCTTTTATCAC", "tT-ACAGcTGCaTTTATCAC", "GG");
```
### Command line (Rust binary)
```bash
# ./cfd_score_calculator <20nt_aligned_spacer_sequence> <20nt_aligned_protospacer_sequence> <last_2nt_of_PAM>
./cfd_score_calculator CTAACAGTTGCTTTTATCAC tT-ACAGcTGCaTTTATCAC GG
```

### Python
```python
# Import:
from cfd_score_calculator import calculate_cfd
# Usage:
# calculate_cfd(spacer, protospacer, pam)
# Example:
calculate_cfd("CTAACAGTTGCTTTTATCAC", "tT-ACAGcTGCaTTTATCAC", "GG")
```
### Command line (Python)
```bash
# Usage:
# python3 cfd_score_calculator.py <20nt_aligned_spacer_sequence> <20nt_aligned_protospacer_sequence> <last_2nt_of_PAM>
# Example:
python3 cfd_score_calculator.py CTAACAGTTGCTTTTATCAC tT-ACAGcTGCaTTTATCAC GG
```

## References
- CFD score was published in [Doench et al (2016) *Nature Biotechnology*](https://doi.org/10.1038/nbt.3437) and the scoring matrices & Python code are adapted from this original publication
- The mismatch scoring matrix was modified to allow for DNA/RNA bulges (intermediate gaps in the alignment) based on the empirical data in Table S19
    - As implemented in [CRISPRme](https://github.com/pinellolab/CRISPRme/blob/main/PostProcess/mismatch_score.pkl) ([Cancellieri, Zeng, Lin et al. (2023) *Nature Genetics*](https://doi.org/10.1038/s41588-022-01257-y))
- To allow for a gap in the alignment at the most PAM-distal nucleotide (for which no empirical data is available), this case is supported with no penalty (a conservative choice to maximize sensitivity)
