
import sys
Protein = sys.argv[3]
pdb_ouput = sys.argv[4]

cmd.load(Protein, 'prot')
cmd.create('prot2', 'polymer')

cmd.save(pdb_ouput, 'prot2', state=0)
cmd.quit()
