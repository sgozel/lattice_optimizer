"""
Copyright 2026 Samuel GOZEL, GNU GPLv3
"""

import numpy as np

BRAVAIS_VECS = {
    'triangle': {'a1': np.array([1.0, 0.0]), 
                 'a2': np.array([0.5, np.sqrt(3)/2])},
    'square': {'a1': np.array([1.0, 0.0]), 
               'a2': np.array([0.0, 1.0])}
    }

SIMU_TORUS = {
    'triangle': {
        int(8): {
            int(10208): {'t1': np.array([1.0, 2.0]),
                         't2': np.array([3.0, -2.0])}},
        int(12): {
            int(20206): {'t1': np.array([2.0, 2.0]),
                         't2': np.array([2.0, -4.0])}},
        int(14): {
            int(10314): {'t1': np.array([1.0, 3.0]),
                         't2': np.array([4.0, -2.0])}},
        int(16): {
            int(40004): {'t1': np.array([4.0, 0.0]),
                         't2': np.array([0.0, 4.0])}},
        int(20): {
            int(10520): {'t1': np.array([1.0, 5.0]),
                         't2': np.array([4.0, 0.0])},
            int(10820): {'t1': np.array([2.0, -4.0]),
                         't2': np.array([3.0, 4.0])}},
        int(21): {
            int(10421): {'t1': np.array([1.0, 4.0]),
                         't2': np.array([5.0, -1.0])},
            int(10821): {'t1': np.array([2.0, -5.0]),
                         't2': np.array([3.0, 3.0])}},
        int(24): {
            int(10424): {'t1': np.array([1.0, 4.0]),
                         't2': np.array([5.0, -4.0])}},
        int(28): {
            int(20414): {'t1': np.array([2.0, 4.0]),
                         't2': np.array([6.0, -2.0])}}
        },
    'square': {
        int(25): {
            int(50005): {'t1': np.array([5.0, 0.0]),
                         't2': np.array([0.0, 5.0])}}
        }
    }
