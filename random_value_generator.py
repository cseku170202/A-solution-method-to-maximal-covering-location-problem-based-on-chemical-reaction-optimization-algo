import matplotlib.pyplot as plt
import numpy

data = []

def start():
    Random()


def Random():
    while True:
        standard_deviation = 15
        mean = 80
        standard_deviation = float(standard_deviation)
        mean = float(mean)
        size = 700 #number of demands to be generated
        size = int(size)
        new_array = numpy.random.uniform(mean, standard_deviation, size)
        print("Generated value for x", new_array)

        fh = open('B700_demand.txt', 'w')
        numpy.savetxt('B700_demand.txt', new_array)
        plt.hist(new_array, 10)
        plt.show()
        fh.close()
        break

Random()
