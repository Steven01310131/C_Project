import matplotlib as plt
import matplotlib.pyplot as plt
import numpy as np


# f = open('Plot.txt', 'r')

x = [1,2,3]
y = [45.48294,9.47978,6.31864]
# for line in f:
#     x.append(float(line.split(',')[0]))
#     y.append(float(line.split(',')[1]))
plt.plot(x, y)
plt.xlabel('steps')
plt.ylabel('Time in seconds')
# plt.title('Plot of time vs number of points')
# fit a square line to the data
# z = np.polyfit(x, y, 2)
# p = np.poly1d(z)
# plt.plot(x, p(x), 'r--')
plt.show()