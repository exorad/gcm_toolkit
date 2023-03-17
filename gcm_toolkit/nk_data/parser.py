from glob import glob

files = glob('*.dat')

out = open('out.txt', 'w')

for fil in files:
    data = open(fil, 'r')
    out.write('    "' + fil[:-8] + '": [\n')
    for line_tot in data.readlines()[5:]:
        line = line_tot.split()
        out.write('        [' + str(line[0]) + ', ' + str(line[0]) + ', '  + str(line[0]) + '],\n')
    out.write('        ],\n\n')
