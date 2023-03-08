def fibonacci_number(n: int):
    prev = curr = 1
    n -= 2

    while (n > 0):
        prev, curr = curr, prev + curr
        n -= 1

    return curr


print(fibonacci_number(3))