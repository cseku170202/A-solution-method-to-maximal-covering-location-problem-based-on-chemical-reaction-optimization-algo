import random
import numpy as np

list = []
a = 0
arr = np.array([])

while a<=700:
      a = a + 1
      b = random.randint(0, 30)
      file = open('test.txt', 'w')
      print(arr)

arr = str(arr)
for i in arr:
  file.write(i)
file.close()




