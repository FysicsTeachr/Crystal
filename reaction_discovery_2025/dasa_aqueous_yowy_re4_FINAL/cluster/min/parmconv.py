import parmed as pmd

structure = pmd.load_file('optimized.pdb')
structure.save('optimized.rst7', format='rst7')
