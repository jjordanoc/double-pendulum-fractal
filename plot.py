import fortranformat as ff
import numpy as np
import matplotlib.pyplot as plt

reader =ff.FortranRecordReader('(3g24.16)')

data = []

with open("data.txt", "r") as file:
    for line in file.readlines():
        data.append(reader.read(line))

d = np.array(data, dtype=np.float64)


plt.imshow(d[:,2].reshape((25, 25)), interpolation="quadric")
# plt.plot(d[:,0], d[:,1])
plt.show()