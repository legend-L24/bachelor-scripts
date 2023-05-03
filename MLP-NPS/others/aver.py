import numpy as np
>>> with open("./rdf.txt") as f:
...     lines = (line for line in f if len(line.split())>3)
...     arr = np.genfromtxt(lines)

>>> for i in range(1, 155):
...     arr_0 = arr_0 + arr[i*500: i*500+500, :]
