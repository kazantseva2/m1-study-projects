import matplotlib.pyplot as plt

# Чтение данных из файла
x_values = []
y_values = []
x_values_local = []
y_values_local = []

with open('/home/varvara/Учеба/magistr_1sem/ПОД/data.txt', 'r') as file:
    for line in file:
        x, y = map(float, line.split())
        x_values.append(x)
        y_values.append(y)
        if (x <= 10):
            x_values_local.append(x)
            y_values_local.append(y)
        

# Первый рисунок
plt.figure(1)  # Создание первой фигуры
plt.plot(x_values, y_values, marker='o', linestyle='-', color='b')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('График зависимости Y от X (полный)')
plt.grid(True)
plt.show()

# Второй рисунок
plt.figure(2)  # Создание второй фигуры
plt.plot(x_values_local, y_values_local, marker='o', linestyle='-', color='r')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('График зависимости Y от X (часть для x <= 2)')
plt.grid(True)
plt.show()
