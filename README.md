# Cutting frequency determination (CFD) score calculator
- A simple, standalone CFD score calculator that supports gaps in RNA-DNA alignment at any non-PAM position
- Implemented in both Rust (faster) and Python (original)
- CFD score is an empirical measure of the likelihood of a given guide RNA (spacer) sequence complexed with SpCas9 cleaving a given DNA (protospacer) sequence. CFD scores are helpful for prioritizing candidate off-target sites for further evaluation. They range from 0 (minimal cutting potential) to 1 (for on-target site)

## Input
- 20 nt aligned spacer (gRNA) sequence
- 20 nt aligned protospacer (DNA) sequence
- Last 2 nt of PAM sequence

(Case insensitive with gaps indicated by `-`)

## Usage
### Rust
```rust
// Import:
mod cfd_score_calculator;
// Usage:
// cfd_score_calculator::calculate_cfd(<20nt_aligned_spacer_sequence>, <20nt_aligned_protospacer_sequence>, <last_2nt_of_PAM>);
// Example:
cfd_score_calculator::calculate_cfd("CTAACAGTTGCTTTTATCAC", "tT-ACAGcTGCaTTTATCAC", "GG");
```
### Command line (Rust binary)
```bash
# Usage:
chmod +x cfd_score_calculator
# ./cfd_score_calculator <20nt_aligned_spacer_sequence> <20nt_aligned_protospacer_sequence> <last_2nt_of_PAM>
# Example:
./cfd_score_calculator CTAACAGTTGCTTTTATCAC tT-ACAGcTGCaTTTATCAC GG
```

### Python
```python
# Import:
from cfd_score_calculator import calculate_cfd
# Usage:
# calculate_cfd(<20nt_aligned_spacer_sequence>, <20nt_aligned_protospacer_sequence>, <last_2nt_of_PAM>)
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

## Implementation notes & references
- CFD score was published in [Doench et al (2016) *Nature Biotechnology*](https://doi.org/10.1038/nbt.3437). The scoring matrices and Python code are adapted from the original publication
- The mismatch scoring matrix was modified to allow for DNA/RNA bulges (intermediate gaps in the alignment) based on the empirical data in Table S19 of Doench et al. 2016
    - As implemented in [CRISPRme](https://github.com/pinellolab/CRISPRme/blob/main/PostProcess/mismatch_score.pkl) ([Cancellieri, Zeng, and Lin et al. (2023) *Nature Genetics*](https://doi.org/10.1038/s41588-022-01257-y))
- To allow for CFD score calculation when there is a gap in the alignment at the most PAM-distal nucleotide (for which no empirical data is available), this case is supported here with no penalty (a conservative implementation choice to maximize sensitivity when nominating candidate off-target sites)