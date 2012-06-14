Simulating disease spread in social networks
============================================
Final project of the "Introduction to Natural Language Processing course"
(ETH Zurich, Spring Semester 2012).

Author
------
Alkis Gkotovos (<alkisg@student.ethz.ch>)

Description and usage example
-----------------------------
Python script `diffuse.py` implements the SIMR model for disease spread.
If run as the main script, it generates all plots used in the report.
Alternatively, it can be imported and used in other modules, as shown
in the example below.

```python
# Import SIMR module and matplotlib for plotting
import diffuse
import matplotlib.pyplot as plt

# Create a new Erdos-Renyi graph with n = 50 and p = 0.1
# with 10% initially infected nodes (randomly chosen)
g = diffuse.Simr(gtype='ER', gparam=[50, 0.1], k=0.1)

# Plot graph state
g.plot()
plt.show()

# Simulate for 20 steps and plot again
for i in range(20):
    g.step()
g.plot()
plt.show()
```