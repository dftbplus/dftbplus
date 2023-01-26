import numpy as np

real_ref = np.array([[[1.73421713, 1.68107958, 1.17731844, 2.],
                      [0.70369243, 0.        , 0.        , 0.],
                      [0.70369243, 0.        , 0.        , 0.]]],
                    dtype=np.float64)

integer_ref = np.array([[[1, 2, 3, 4, 5, 6]]], dtype=np.int64)

complex_ref = np.array([[1.1 + 2.2j, 3.3 + 4.4j, 5.5 + 6.6j],
                        [7.7 + 8.8j, 9.9 + 10.1j, 11.11 + 12.12j],
                        [13.13 + 14.14j, 15.15 + 16.16j, 17.17 + 18.18j]],
                       dtype=np.clongdouble)

logical_ref = np.array([[[False]],
                        [[True]],
                        [[True]],
                        [[False]],
                        [[True]],
                        [[False]]],
                       dtype=bool)

scalar_ref = np.array([0.173421713306389E+001], dtype=np.float64)
