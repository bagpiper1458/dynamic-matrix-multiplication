import random
import unittest
from dynamic_matrix.audit4 import make_audit
from dynamic_matrix import gf_mul


class AuditTests(unittest.TestCase):
    def test_formula_mapping_and_counts(self):
        for kind, coefficient_xors, retained, mults in [('toeplitz',11,9,9),('circulant',5,5,9),('unit-circulant',5,5,5)]:
            audit = make_audit(kind)
            self.assertEqual(audit['counts']['coefficient_xors'], coefficient_xors)
            self.assertEqual(audit['counts']['retained_coefficients'], retained)
            self.assertEqual(audit['counts']['multiplications'], mults)
            self.assertEqual(audit['counts']['vector_xors'],15)
            self.assertEqual(audit['counts']['identity_copies'],9-mults)
            self.assertEqual(len(audit['formula_terms']),16)
            self.assertEqual([[t['product'] for t in audit['formula_terms'] if t['output']==f'y{i}'] for i in range(4)],
                             [['p0','p1','p3','p4'],['p0','p2','p3','p5'],['p0','p1','p6','p7'],['p0','p2','p6','p8']])
            self.assertTrue(all(c['verified'] for c in audit['certificate']))
            self.assertEqual(audit['counts']['peak_slots'],9)

    def test_register_replay_without_freed_values(self):
        rng = random.Random(20261008)
        for kind in ['toeplitz','circulant','unit-circulant']:
            for _ in range(50):
                column = [rng.randrange(256) for _ in range(4)]
                if kind=='unit-circulant':
                    column[-1]=1^column[0]^column[1]^column[2]
                row = [column[0],*[rng.randrange(256) for _ in range(3)]]
                vector=[rng.randrange(256) for _ in range(4)]
                audit=make_audit(kind,column,row,vector)
                persistent={c['id']:c['value'] for c in audit['persistent_coefficients']}
                raw={f'b{i}':b for i,b in enumerate(vector)}
                registers=[None]*9
                for step in audit['vector_steps']:
                    def read(name):
                        if name in raw: return raw[name]
                        if name in persistent: return persistent[name]
                        return registers[step['live_before'][name]]
                    operands=[read(dep) for dep in step['deps']]
                    result=operands[0] if step['op']=='copy' else (operands[0]^operands[1] if step['op']=='xor' else gf_mul(*operands))
                    for name in step['released']:
                        registers[step['live_before'][name]]=None
                    registers[step['slot']]=result
                    self.assertEqual(result,step['value'])
                last=audit['vector_steps'][-1]['live_after']
                self.assertEqual([registers[last[f'y{i}']] for i in range(4)],audit['expected'])


if __name__=='__main__': unittest.main()
