from multiprocessing import Pool
import time


COUNT = 99999999
def countdown(n):
    while n>0:
        n -= 1

"""
start = time.time()
countdown(COUNT)
end = time.time()
print("Single thread time taken in seconds:", end - start)
"""

if __name__ == '__main__':
    pool = Pool(processes=3)
    start = time.time()
    r1 = pool.apply_async(countdown, [COUNT//3])
    r2 = pool.apply_async(countdown, [COUNT//3])
    r3 = pool.apply_async(countdown, [COUNT//3])
    pool.close()
    pool.join()
    end = time.time()
    print('Multi thread time taken in seconds:', end - start)

start = time.time()
countdown(COUNT)
end = time.time()
print("Single thread time taken in seconds:", end - start)