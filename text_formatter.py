import csv

data = open('valid_Lua_2600_cmds.txt', 'r')
linez = data.readlines()
new_linez = []
for line in linez:
    for i in range(len(line)):
        if line[-1].isalpha():
            break
        if (line[i] == '.' and line[i+1] == '.'):
            new_linez.append(line[0:i])
            break

datanew = open('Lua_2600_cmds.txt', 'w')
for i in range(len(new_linez)):
    datanew.writelines(new_linez[i])
    datanew.writelines('\n')
