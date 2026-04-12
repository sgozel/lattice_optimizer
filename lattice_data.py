"""
Copyright 2026 Samuel GOZEL, GNU GPLv3
"""

import numpy as np

BRAVAIS_VECS = {
    'triangle': {
        int(16):  {'a1': np.array([1.0, 0.0]), 
                   'a2': np.array([0.5, np.sqrt(3)/2])},
        int(20): {'a1': np.array([1.0, 0.0]), 
                  'a2': np.array([0.5, np.sqrt(3)/2])},
        int(24): {'a1': np.array([1.0, 0.0]), 
                  'a2': np.array([0.5, np.sqrt(3)/2])},
        int(28): {'a1': np.array([1.0, 0.0]),
                  'a2': np.array([0.5, np.sqrt(3)/2])}
        },
    'square': {
        int(25): {'a1': np.array([1.0, 0.0]), 
                  'a2': np.array([0.0, 1.0])}
        }
    }

SIMU_TORUS = {
    'triangle': {
        int(16): {'t1': np.array([4.0, 0.0]),
                  't2': np.array([0.0, 4.0])},
        int(20): {'t1': np.array([2.0, -4.0]),
                  't2': np.array([3.0, 4.0])},
        int(24): {'t1': np.array([1.0, 4.0]),
                  't2': np.array([5.0, -4.0])},
        int(28): {'t1': np.array([2.0, 4.0]),
                  't2': np.array([6.0, -2.0])}
        },
    'square': {
        int(25): {'t1': np.array([5.0, 0.0]),
                  't2': np.array([0.0, 5.0])}
        }
    }
