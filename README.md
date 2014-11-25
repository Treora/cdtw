Constrained Dynamic Time Warping
================================

This is a small Python module, written in C, implementing the cDTW similarity measure between two sequences of numbers. This implementation is specifically optimised for using the Sakoe-Chiba band as the constraint on the time shift, and provides just one function:

`cdtw_sakoe_chiba(sequence1, sequence2, r)`

The function returns a float representing the dissimilarity between the sequences. The sequences must be one-dimensional NumPy arrays containing float64 values. This implementation requires the sequences to be of equal length, although in theory this could be alleviated slightly, allowing a difference of `r`.

The width of the band is defined by the third argument `r`, which sets the maximum number of shifted time-steps. Oftentimes, `r` is specified as a percentage of the sequence length, but the conversion is left to you. Passing `0` for `r` will return the L1 norm (taxicab/Manhattan distance), whereas passing `len(sequence)` would result in doing unconstrained DTW.

Note that no normalisation is performed on the sequences, which may be desired.

### Build and run

For convenience, running `make` builds the module, installs it in the current directory, and runs the test example found below. Run `python setup.py install` to build the module and install it system-wide.

Due to [issue #2](https://github.com/Treora/cdtw/issues/2), you need to have NumPy installed before building this module (e.g. via `pip install numpy`), this is not done automatically.

### Example usage

```python
import numpy as np
from cdtw import cdtw_sakoe_chiba

a = np.array([1,2,2,2,3,4,5,6,7], dtype="float64")
b = np.array([1,2,3,4,5,6,6,6,7], dtype="float64")
d = [cdtw_sakoe_chiba(a, b, r) for r in [0,1,2]]
print d
```

Running this code prints `[8.0, 4.0, 0.0]`. See also the example in `cdtw_example.py`.

### Metastuff

Please feel free to test and improve the code, add more functions to the module (other constraints, lower bounds), and report any bugs you may encounter. There seem to be plenty of other implementations of DTW available; this one was mostly written for the fun of it, and to create a fast version optimised for using the Sakoe-Chiba band constraint.

This is free and unencumbered software released into the public domain.
