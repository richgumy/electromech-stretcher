data = open('Lua_cmds_good.txt', 'r')

linez = data.readlines()

datanew = open('lua_python_code.txt', 'w')
lines_split = []
for line in linez:
    line = line.strip('\n* ')
    lines_split.append(line.split('.'))

# print(lines_split)

same_class = 0
for i in range(len(lines_split)): # [[beeper, beep()],[beeper,enable],[bit,bitand() ]]
    print(i)
    if lines_split[i][0] == lines_split[i-1][0]:
        same_class = 1
    else:
        same_class = 0
    for j in range(len(lines_split[i])):
        print(j)
        for tabs in range(j):
            datanew.writelines('\t') # writes number of tabs as the index in line list
        if (j == (len(lines_split[i])-1)) and (lines_split[i][j][-1] == ')'):
            # If its a function (@ EOL) write code in format:
            # def clear(self):
            #     self.connection.write("display.clear()")
            datanew.writelines('def ')
            datanew.writelines(lines_split[i][j].strip('()'))
            datanew.writelines('(self):\n')
            for tabs in range(j):
                datanew.writelines('\t') # writes number of tabs as the index in line list
            datanew.writelines('\tself.connection.write("')
            for k in (range(len(lines_split[i])-1)):
                datanew.writelines(lines_split[i][k])
                datanew.writelines('.')
            datanew.writelines(lines_split[i][j])
            datanew.writelines('()")\n') # can got through and manually fill in function parameters

        if (j == (len(lines_split[i])-1)) and (lines_split[i][j][-1] != ')'):
            # If its a val (@ EOL) write code in format:
            # def func(self, value=None):
            #     if value != None:
            #         self.connection.write("display.smua.measure.func = %d" % value)
            #     else:
            #         self.connection.write("meas_func = display.smua.measure.func")
            #         return self.connection.query("print(meas_func)")
            datanew.writelines('def ')
            datanew.writelines(lines_split[i][j])
            datanew.writelines('(self, value=None):\n')

            for tabs in range(j):
                datanew.writelines('\t') # writes number of tabs as the index in line list
            datanew.writelines('\tif value != None:\n')

            for tabs in range(j):
                datanew.writelines('\t') # writes number of tabs as the index in line list
            datanew.writelines('\t\tself.connection.write("')
            for k in (range(len(lines_split[i])-1)):
                datanew.writelines(lines_split[i][k])
                datanew.writelines('.')
            datanew.writelines(lines_split[i][j])
            datanew.writelines(' = %d" % value)\n') # can got through and manually fill in function parameters

            for tabs in range(j):
                datanew.writelines('\t') # writes number of tabs as the index in line list
            datanew.writelines('\telse:\n')

            for tabs in range(j):
                datanew.writelines('\t') # writes number of tabs as the index in line list
            datanew.writelines('\t\tself.connection.write("')
            for k in (range(len(lines_split[i]))):
                datanew.writelines(lines_split[i][k])
            datanew.writelines(' = ')
            for k in (range(len(lines_split[i])-1)):
                datanew.writelines(lines_split[i][k])
                datanew.writelines('.')
            datanew.writelines(lines_split[i][j])
            datanew.writelines('")\n')

            for tabs in range(j):
                datanew.writelines('\t') # writes number of tabs as the index in line list
            datanew.writelines('\t\treturn self.connection.query("print(')
            for k in (range(len(lines_split[i]))):
                datanew.writelines(lines_split[i][k])
            datanew.writelines(')")\n')

        else:
            # if it's an obj we want the format:
            # class display:
            #     class smua:
            #         class measure:
            #           ...
            if same_class == 0:
                for tabs in range(j):
                    datanew.writelines('\t') # writes number of tabs as the index in line list
                datanew.writelines('class ')
                datanew.writelines(lines_split[i][j])
                datanew.writelines(':\n')



# example:
# class display:
#     def clear(self):
#         self.connection.write("display.clear()")
#     def settext(self, text: str):
#         self.connection.write("display.settext('%s')" % text)
#         print("display.settext('%s')" % text)
#     class smua:
#         class measure:
#             def func(self, value=None):
#                 if value != None:
#                     self.connection.write("display.smua.measure.func = %d" % value)
#                 else:
#                     self.connection.write("meas_func = display.smua.measure.func")
#                     return float(self.connection.query("print(meas_func)"))
