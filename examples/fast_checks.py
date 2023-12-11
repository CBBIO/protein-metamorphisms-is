from Bio import ExPASy, SwissProt

handle = ExPASy.get_sprot_raw("A0A1S4F020")
record = SwissProt.read(handle)

for cr in record.cross_references:
    if cr[0]=="PDB":
        print(cr)


a = "1/2/0/4/5/3/7/A/6/C/D/B/E/F/G/I/J/H/L/M/K/O/P/N/R/Q/T/S/V/U/X/W/Z/Y/b/a/d/c/f/e/h/g/j/i/l/k/n/m/p/o/r/q/t/s/v/u/x/w/z/y"
a = set(a.split('/'))

b = "0/1/2/3/4/5/6/7/A/B/C/D/E/F/G/H/I/J/K/L/M/N/O/P/Q/R/S/T/U/V"
b = set(b.split('/'))

print(a)
print(b)
