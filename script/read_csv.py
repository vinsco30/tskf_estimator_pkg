
# import csv

# with open('/home/vinsco/Desktop/test2.csv', 'r', encoding='utf-16') as csv_file:
#     csv_reader = csv.reader(csv_file)
#     line_cnt = 0
#     idx = 0
#     for row in csv_reader:
#         if line_cnt ==0:
#             print(f'Column names are {", ".join(row)}')
#             line_cnt += 1
#         else:
#             print(f'Fine simulazione, il numero di righe della simulazione è: {line_cnt}')
#             else:
#                 line_cnt += 1
#     print(f'Processed {line_cnt} lines.')

import csv
with open('/home/vinsco/Desktop/test2.csv', 'r', encoding='utf-16') as file:
    reader = csv.reader(file, delimiter = ',')
    r=0
    c = 0
    for r=row in reader:
        
#     # print(reader)
