

with open('output.txt', 'r+') as file:
    line = file.readline()
    line_num = 0
    while line:
        buffer = line
        line = file.readline()
        # line_num += 1
        # file.seek(line_num)
        file.write(buffer)

    # file.write(data)

print("Данные успешно записаны в файл output.txt")