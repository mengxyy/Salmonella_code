import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import t

# Load gene matrix and remove first column
df = pd.read_excel('input.xlsx')
df_genes = df.iloc[:, 1:]

# Sampling configurations
sample_sizes = [10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000, 200000, 300000, 400000, 460000]
repetitions = 5
results = {}

# Random sampling and count unique ARGs (value > 0)
for n in sample_sizes:
    counts = []
    for i in range(repetitions):
        sample = df_genes.sample(n=n, replace=False, random_state=i)
        unique_genes = (sample > 0).any().sum()
        counts.append(unique_genes)
    results[n] = counts


# Fit logarithmic curve to mean values
means = results_df.mean().values
stds = results_df.std().values
x_data = np.array(sample_sizes)

def log_func(x, a, b):
    return a * np.log(x) + b

popt, _ = curve_fit(log_func, x_data, means)
r2 = 1 - np.sum((means - log_func(x_data, *popt))**2) / np.sum((means - means.mean())**2)

# Plot
x_fit = np.linspace(x_data.min(), x_data.max(), 500)
y_fit = log_func(x_fit, *popt)

plt.figure(figsize=(8, 6))
plt.errorbar(x_data, means, yerr=stds, fmt='ko', capsize=5, label='Observed')
plt.plot(x_fit, y_fit, 'r-', label=f'Fit: $y = {popt[0]:.2f} \\ln(x) + {popt[1]:.2f}$, $R^2 = {r2:.3f}$')
plt.xlabel('Number of Strains')
plt.ylabel('Number of ARG Types')
plt.legend()
plt.grid(True)
plt.show()
