delta = 1e-8
eps = 1e-5

def objective(func, x):
    return eval(func, {"x": x})


def dichotomy(begin: float, end: float, func: str, eps: float, delta: float):
    """
    Finds an interval that contains the minimum of a function
    Also creates a file with csv table of the following format:

    +------+------+-------+------------+-------------+----------+--------+-------------------+------------------------------------------------+
    | iter | left | right | value_left | value_right | begin(i) | end(i) | end(i) - begin(i) | end(i - 1) - begin(i - 1) /  end(i) - begin(i) |
    +------+------+-------+------------+-------------+----------+--------+-------------------+------------------------------------------------+

    :param begin: begin of initial interval
    :param end: end of initial interval
    :param func: function expression in Python format, must be
                 valid for standard eval() function
    :param eps: precision
    :param delta: value of deflection from the middle of interval

    :return: void
    """

    iters = 0

    while (abs(end - begin) > eps):

        iters += 1

        left = (begin + end - delta) / 2
        right = (begin + end + delta) / 2

        value_left = objective(func, left)
        value_right = objective(func, right)

        if (value_left < value_right):
            end = right
        elif (value_right < value_left):
            begin = left
        else:
            begin = left
            end = right

    print(f"{begin:.8e}, {end:.8e}, {iters}")

def golden_section(begin, end, func):

    golden_number = 0.381966011

    left = begin + golden_number * (end - begin)
    right = end - golden_number * (end - begin)

    value_left = objective(func, left)
    value_right = objective(func, right)

    while (abs(end - begin) > eps):
        if (value_left < value_right):
            end = right
            right = left
            value_right = value_left
            left = begin + golden_number*(end - begin)
            value_left = objective(func, left)
        else:
            begin = left
            left = right
            value_left = value_right
            right = end - golden_number*(end - begin)
            value_right = objective(func, right)

    print(f"{begin:.8e}, {end:.8e}")


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