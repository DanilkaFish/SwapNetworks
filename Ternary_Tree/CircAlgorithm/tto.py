a = [92, 94, 86, 97, 95, 98, 97, 97, 88, 84, 88, 100, 92, 84, 86, 80, 74, 78, 82, 84, 94, 88, 88, 98]
b = {i: a.count(i) for i in a}
print(b)
import matplotlib.pyplot as plt
import numpy as np


plt.hist(a, align='left', rwidth=0.9)  # density=False would make counts
plt.ylabel('Probability')
plt.xlabel('Data');
plt.savefig('histogram.png')