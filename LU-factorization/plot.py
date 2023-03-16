import matplotlib as plt
import matplotlib.pyplot as plt
import numpy as np


f = open('pivoting.txt', 'r')
# f1 = open('Plot1.txt', 'r')
# f2= open('Plot2.txt', 'r')

x = []
y = []
for line in f:
    x.append(float(line.split(',')[1]))
    y.append(float(line.split(',')[0]))
print(x)
# x1= []
# y1= []
# for line in f1:
#     x1.append(float(line.split(',')[1]))
#     y1.append(float(line.split(',')[0]))
# x2= []
# y2= []
# for line in f2:
#     x2.append(float(line.split(',')[1]))
#     y2.append(float(line.split(',')[0]))
plt.plot(x, y)
# plt.plot(x1, y1,label="Flags")
# plt.plot(x2, y2,label="Final version")
plt.xlabel('Size of NxN matrix')
plt.ylabel('Time in seconds')
# fit a square line to the data
# z = np.polyfit(x, y, 3)
# p = np.poly1d(z)
# plt.plot(x, p(x), 'r--',label="theoretical")
plt.legend()
plt.show()