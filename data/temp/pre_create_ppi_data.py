import csv

f = open('test.csv', 'rU')
csv_f = csv.reader(f)
csv_list = list(csv_f)

res = []
for i in csv_list:
    tmp = (i[0], i[1])
    if ((i[1], i[0]) not in res)  and (tmp not in res):
        res.append(tmp)

result = open('data/ppi1.py','w')
result.write('ppi = \n')
result.write(str(res))
