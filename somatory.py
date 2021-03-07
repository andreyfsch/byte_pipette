n = 348

def calc(n):
    sum = 1
    for i in range(n, 1, -1):
        print(i)
        sum += i*(i-1)
    return sum

print(calc(n))