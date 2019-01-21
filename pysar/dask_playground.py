import dask
import math

from dask.distributed import Client, as_completed
client = Client(processes=False)

d = {}
e = {}
futures = []

def f(i, j):
    return (i, j), math.pow(i ** 50 + j ** 50, 1 / 50)

for i in range(100):
    for j in range(100):
        future = client.submit(f,i, j)
        futures.append(future)

for future in as_completed(futures):
    (i, j), res = future.result()
    d[(i, j)] = res

print(d)
