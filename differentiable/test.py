#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np


X = np.load("data/x_data.npy")
Y = np.load("data/y_data.npy")
g = np.load("data/g_data.npy")
derivatives_data = np.load("data/dgdk_data.npy.npz")
dgdk_values, dgdk_bend_values, dgdk_values_finite, dgdk_bend_values_finite = [derivatives_data[dfile] for dfile in derivatives_data.files]

magnitudes = np.sqrt(dgdk_bend_values**2 + dgdk_values**2)

plt.contourf(X, Y, g, levels=300)
plt.quiver(X, Y, dgdk_values/magnitudes, dgdk_bend_values/magnitudes, color="blue")
plt.quiver(X, Y, dgdk_values_finite/magnitudes, dgdk_bend_values_finite/magnitudes , color="red")
plt.show()
