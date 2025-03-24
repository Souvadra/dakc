import numpy as np
from math import log2, ceil

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, show
import seaborn as sns
from matplotlib.pyplot import plot, show

import math

# Set the ggplot style using Seaborn
sns.set_theme(style="darkgrid", palette="dark")

# Initialize the colors array using Seaborn's color_palette function
colors = sns.color_palette("dark")

# Set the font properties globally to use Helvetica
# plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 10  # Example: 10, 12, 14
plt.rcParams['font.weight'] = 'bold'  # Example: 'normal', 'bold'
