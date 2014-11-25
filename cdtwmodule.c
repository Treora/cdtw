
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>


#define E_OUT_OF_MEMORY (-1)

#define validate_input(assertion, message) \
{ \
    if (!(assertion)) { \
        PyErr_SetString(PyExc_ValueError, message); \
        return NULL; \
    } \
}

#ifdef max
#undef max
#endif
#define max(a, b) (a>b ? a : b)

#ifdef min
#undef min
#endif
#define min(a, b) (a<b ? a : b)

#define min3(a,b,c) (min(min(a,b),c))

double
_l1_norm(PyArrayObject* arr1, PyArrayObject* arr2, int seq_length)
{
    double x1, x2;
    double distance = 0;
    int i;
    for (i=0; i<seq_length; i++) {
        x1 = * (double*) PyArray_GETPTR1(arr1, i);
        x2 = * (double*) PyArray_GETPTR1(arr2, i);
        distance += fabs(x1-x2);
    }
    return distance;
}


double
_cdtw_sakoe_chiba(PyArrayObject* arr1, PyArrayObject* arr2, int seq_length, int r)
{
    double x1, x2;
    int seq1_length, seq2_length;
    // For now, both sequences have to be of the same length. For low values of r this restriction would be implied anyway.
    seq1_length = seq2_length = seq_length;

    // Note: To save a little memory, our 'column' is not a full column of the DTW matrix, but only the part within r steps from the diagonal (although if 2*r+1 > seq2_length, this approach actually costs more memory). This makes the formulas a bit different than usual. For example, curr_column[1+r] always represents the element on the diagonal.
    // Columns need r+1+r elements to cover the range of the allowed time skew on both sides of the diagonal of the DTW matrix. Two extra elements are fixed to infinity to simplify finding best_value further below.
    int column_size = 1+r+1+r+1;
    double (*prev_column)[column_size] = malloc(column_size*sizeof(double));
    if (prev_column == NULL) {
        return E_OUT_OF_MEMORY;
    }
    double (*curr_column)[column_size] = malloc(column_size*sizeof(double));
    if (curr_column == NULL) {
        return E_OUT_OF_MEMORY;
    }

    // In the DTW matrix the items left of the first column must be infinite, as must values below the bottom row.
    int i;
    for (i=0; i<column_size; i++) {
        (*prev_column)[i] = INFINITY;
        (*curr_column)[i] = INFINITY;
    }
    // We do of course need a place to start from. This forces the first matrix element to become |arr1[0]-arr2[0]|
    (*prev_column)[r+1] = 0;

    int top_valid_row=0, bottom_valid_row;
    int col, row;
    for (col=0; col<seq1_length; col++) {
        // Pick range of rows, ensuring 0 <= j < seq2_length
        bottom_valid_row = max(1, -col+r+1);
        top_valid_row = min(column_size-2, seq2_length-1-col+r+1);

        // Fill up this column
        for (row=bottom_valid_row; row<=top_valid_row; row++) {
            // Compute mapping to the index in the full DTW matrix.
            int i = col;
            int j = col - r + row-1;
            // Read values from the numpy array
            x1 = * (double*) PyArray_GETPTR1(arr1, i);
            x2 = * (double*) PyArray_GETPTR1(arr2, j);

            // Compute their absolute error
            double error = fabs(x1-x2);

            // Determine via which path to get here with lowest cumulative error
            double left_value = (*prev_column)[row+1];
            double left_under_value = (*prev_column)[row];
            double under_value = (*curr_column)[row-1];
            double best_value = min3(left_under_value, left_value, under_value);
            // Remember the best possible cumulative error for this point
            (*curr_column)[row] = error + best_value;
        }
        // Move to the next column, forget the previous one. (*curr_column is only written to, so it need not be cleaned)
        double (*swap_temp)[column_size] = prev_column;
        prev_column = curr_column;
        curr_column = swap_temp;
    }
    // The result can be found in the upper right corner of the DTW matrix
    double result = (*prev_column)[top_valid_row];

    free(prev_column);
    free(curr_column);

    return result;
}


static PyObject*
cdtw_sakoe_chiba(PyObject* self, PyObject* args)
{
    // Inputs
    PyArrayObject* arr1 = NULL;
    PyArrayObject* arr2 = NULL;
    int r;
    // Output
    double dissimilarity;

    int seq_length;

    // Read the arguments from Python into the C variables.
    if (!PyArg_ParseTuple(args, "O!O!i", &PyArray_Type, &arr1, &PyArray_Type, &arr2, &r))
        return NULL;

    // We can operate only on one-dimensional arrays (sequences).
    validate_input(PyArray_NDIM(arr1) == 1, "First sequence is not one-dimensional.");
    validate_input(PyArray_NDIM(arr2) == 1, "Second sequence is not one-dimensional.");
    // Check that the array elements are aligned in memory and of type double. We could convert them ourselves, but leave it to the programmer for now.
    validate_input(PyArray_ISALIGNED(arr1) && PyArray_ISALIGNED(arr2),
        "NumPy arrays should be properly aligned in memory, somehow this is not the case.."
    );
    // We only accept arrays of doubles (dtype='float64') for now.
    validate_input(PyArray_TYPE(arr1)==NPY_DOUBLE, "First sequence is not of type double/float64");
    validate_input(PyArray_TYPE(arr2)==NPY_DOUBLE, "Second sequence is not of type double/float64");
    // We require both arrays to be of the same length for now.
    validate_input(PyArray_DIM(arr1, 0) == PyArray_DIM(arr2, 0),
        "The sequences should be of equal length."
    );

    seq_length = PyArray_DIM(arr1, 0);
    validate_input(r <= seq_length, "The constraint on the time-shift, r, cannot be larger than the sequences' lengths");

    if (r == 0) {
        // For the case r=0, we can easily optimise the algorithm as it is just the L1 norm (sum of absolute errors).
        dissimilarity = _l1_norm(arr1, arr2, seq_length);
    }
    else {
        dissimilarity = _cdtw_sakoe_chiba(arr1, arr2, seq_length, r);
        if (dissimilarity == E_OUT_OF_MEMORY)
            return PyErr_NoMemory();
    }
    return Py_BuildValue("d", dissimilarity);
}

static PyMethodDef cdtwMethods[] = {
    {"cdtw_sakoe_chiba", cdtw_sakoe_chiba, METH_VARARGS,
        "cdtw_sakoe_chiba(sequence1, sequence2, r) -> dissimilarity" "\n"
        "Use constrained Dynamic Time Warping to calculate the dissimilarity between two sequences. Pass it two one-dimensional NumPy arrays of float64 (double) values.\n"
        "\n"
        "The value r determines the size of the time-shift-constraining Sakoe-Chiba band (r=0 gives the L1 norm, r=len(sequence) gives unconstrained DTW)\n"
    },
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initcdtw(void)
{
    (void) Py_InitModule("cdtw", cdtwMethods);
    import_array();
}
