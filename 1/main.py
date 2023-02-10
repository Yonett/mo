delta = 1e-6
eps = 1e-5

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
    print(f"{x1:.8e}, {x2:.8e}")

def golden_section(a, b, func):
    x1 = a + 0.381966011*(b-a)
    x2 = b - 0.381966011*(b-a)
    value1 = objective(func, x1)
    value2 = objective(func, x2)
    while (abs(b-a) > eps):
        if (value1 < value2):
            b = x2
            x2 = x1
            value2 = value1
            x1 = a + 0.381966011*(b-a)
            value1 = objective(func, x1)
        else:
            a = x1
            x1 = x2
            value1 = value2
            x2 = b - 0.381966011*(b-a)
            value2 = objective(func, x2)
    print(f"{x1:.8e}, {x2:.8e}")


def main():
    file = open("funcs.txt")
    for func in file:
        print(func, "precision: ", eps)
        print("Dichotomy:")
        dichotomy(-2, 20, func)
        print("Golden Section:")
        golden_section(-2, 20, func)
        print("")

    file.close()


if __name__ == "__main__":
    main()