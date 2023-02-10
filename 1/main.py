delta = 1e-11
eps = 1e-10

def func(x):
    return (x-1)**2

def objective(func, x):
    return eval(func, {"x": x})

def dichotomy(a, b, func):
    while (abs(b-a) > eps):
        x1 = (a + b - delta) / 2
        x2 = (a + b + delta) / 2
        if (objective(func, x1) < objective(func, x2)):
            b = x2
        elif (objective(func, x2) < objective(func, x1)):
            a = x1
        else:
            a = x1
            b = x2
    print(f"{x1:.10e}, {x2:.10e}")


def main():
    file = open("funcs.txt")
    print("Dichotomy:")
    for func in file:
        dichotomy(-2, 20, func)

    file.close()


if __name__ == "__main__":
    main()