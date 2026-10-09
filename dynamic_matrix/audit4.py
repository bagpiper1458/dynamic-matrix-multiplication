"""Inspectable 4x4 PARAMS reconstruction, with symbolic and numeric certificates.

This module documents one explicit schedule. It does not claim to recover the
lost code's storage layout, or to minimize the number of XORs or live registers.
"""
from collections import Counter
from functools import reduce
from operator import xor

from .core import gf_mul
from .iterative import make_iterative_layout


def symbol(offset):
    return f"a[{offset}]"


def expression(terms):
    return " XOR ".join(sorted(terms)) or "0"


def _tensor(coefficient, vector):
    return {(a, b) for a in coefficient for b in vector}


def _xor_sets(*sets):
    result = set()
    for terms in sets:
        result.symmetric_difference_update(terms)
    return result


def _lifetime(steps, inputs, persistent):
    """Topological use counts; reuse a dead slot only after reading all operands.

    Values held in input/constant storage are excluded. Output slots stay live.
    peak_slots excludes one scalar result register used before slot reuse.
    """
    uses = Counter(dep for step in steps for dep in step['deps'])
    last = {dep: index for index, step in enumerate(steps, 1) for dep in step['deps']}
    outputs = {f'y{i}' for i in range(4)}
    slots, free, next_slot, peak = {}, [], 0, 0
    result = []
    for index, step in enumerate(steps, 1):
        before = dict(slots)
        released = []
        for dep in step['deps']:
            uses[dep] -= 1
            if uses[dep] == 0 and dep in slots and dep not in outputs:
                free.append(slots.pop(dep))
                released.append(dep)
        if free:
            slot = min(free)
            free.remove(slot)
        else:
            slot, next_slot = next_slot, next_slot + 1
        slots[step['id']] = slot
        peak = max(peak, len(slots))
        result.append(dict(step= index, **step, slot=slot, released=released,
                           live_before=before, live_after=dict(slots),
                           last_use='保留輸出' if step['id'] in outputs else last.get(step['id'])))
    if set(slots) != outputs:
        raise AssertionError('liveness schedule did not finish with four outputs')
    if any(dep not in inputs and dep not in persistent and dep not in {s['id'] for s in steps[:i]}
           for i, step in enumerate(steps) for dep in step['deps']):
        raise AssertionError('schedule dependency not available')
    return result, peak


def make_audit(kind='toeplitz', column=None, row=None, vector=None):
    if kind not in ('toeplitz', 'circulant', 'unit-circulant'):
        raise ValueError('unknown audit case')
    circular = kind != 'toeplitz'
    unit = kind == 'unit-circulant'
    column = list(column) if column is not None else ([2,1,1,3] if unit else [2,7,1,3])
    row = list(row) if row is not None else ([2,9,11,13] if not circular else [column[0], *column[:0:-1]])
    vector = list(vector) if vector is not None else [0xDB,0x13,0x53,0x45]
    if len(column) != 4 or len(row) != 4 or len(vector) != 4:
        raise ValueError('audit requires four elements in each input')
    if circular:
        row = [column[0], *column[:0:-1]]
    layout = make_iterative_layout(4, circulant=circular)
    plan = layout.prepare(column, None if circular else row, unit_row_sum=unit)
    # Validate vector using the public evaluator before constructing trace data.
    iterative_result = plan.matvec(vector)
    input_names = [symbol(i) for i in range(4)] if circular else [symbol(i) for i in range(-3,4)]
    raw_values = column if circular else [*reversed(row[1:]), *column]
    names, forms = input_names + ['zero'], [{name} for name in input_names] + [set()]
    values = dict(zip(input_names, raw_values))
    values['zero'] = 0
    coeff_steps = []
    for index, (left, right) in enumerate(layout.coefficient_ops):
        name = f'c{index}'
        expanded = forms[left] ^ forms[right]
        values[name] = values[names[left]] ^ values[names[right]]
        coeff_steps.append(dict(id=name, op='xor', deps=[names[left],names[right]],
                                expression=f'{names[left]} XOR {names[right]}', expanded=expression(expanded),
                                value=values[name]))
        names.append(name)
        forms.append(expanded)
    coefficient_names = [names[i] for i in layout.coefficient_ids]
    coefficient_forms = [forms[i] for i in layout.coefficient_ids]
    persistent = set(coefficient_names)
    for index, step in enumerate(coeff_steps, 1):
        users = [s['id'] for s in coeff_steps if step['id'] in s['deps']]
        step.update(step=index, users=users,
                    retained=step['id'] in persistent,
                    last_use='矩陣不變時保留' if step['id'] in persistent else max(
                        (j for j,s in enumerate(coeff_steps,1) if step['id'] in s['deps']), default=index))

    for i,b in enumerate(vector):
        values[f'b{i}'] = b
    linear = {f'b{i}': {f'b{i}'} for i in range(4)}
    steps = []
    def add(name, left, right):
        values[name] = values[left] ^ values[right]
        steps.append(dict(id=name, op='xor', deps=[left,right], expression=f'{left} XOR {right}', value=values[name]))
        if left in linear and right in linear:
            linear[name] = linear[left] ^ linear[right]
    add('v01','b0','b1')
    add('v23','b2','b3')
    add('vAll','v01','v23')
    add('v13','b1','b3')
    add('v02','b0','b2')
    data_names = ['vAll','v13','v02','v23','b3','b2','v01','b1','b0']
    products, tensors = [], {}
    for i,(a,b) in enumerate(zip(coefficient_names,data_names)):
        name = f'p{i}'
        identity = unit and i in layout.unit_leaves
        values[name] = values[b] if identity else gf_mul(values[a],values[b])
        tensors[name] = _tensor(coefficient_forms[i],linear[b])
        step = dict(id=name, op='copy' if identity else 'multiply', deps=[b] if identity else [a,b],
                    expression=f'{b}（列 XOR = 1）' if identity else f'{a} * {b}', value=values[name])
        steps.append(step)
        products.append(dict(id=name, coefficient=a, coefficient_expanded=expression(coefficient_forms[i]),
                             vector=b, vector_expanded=expression(linear[b]), value=values[name],
                             operation=step['op'], outputs=[]))
    for name,left,right in [('f0','p0','p1'),('f1','p0','p2'),('g0','p3','p4'),('g1','p3','p5'),
                            ('h0','p6','p7'),('h1','p6','p8'),('y0','f0','g0'),('y1','f1','g1'),
                            ('y2','f0','h0'),('y3','f1','h1')]:
        add(name,left,right)
        tensors[name] = tensors[left] ^ tensors[right]
    term_map = []
    for i in range(4):
        for mask in range(4):
            p = i & mask
            offsets = [p-s for s in range(4) if s & mask == s]
            a_set = {symbol(offset % 4 if circular else offset) for offset in offsets}
            b_set = {f'b{((~i)&mask)|s}' for s in range(4) if s & (3^mask) == s}
            matches = [j for j in range(9) if coefficient_forms[j] == a_set and linear[data_names[j]] == b_set]
            if len(matches) != 1:
                raise AssertionError('literal general-formula term did not map uniquely')
            product = products[matches[0]]
            product['outputs'].append(f'y{i}')
            term_map.append(dict(output=f'y{i}', mask=mask, A=expression(a_set), B=expression(b_set), product=product['id']))
    certificate = []
    for i in range(4):
        expected = {(symbol((i-j)%4 if circular else i-j),f'b{j}') for j in range(4)}
        mapped = _xor_sets(*(tensors[t['product']] for t in term_map if t['output']==f'y{i}'))
        if tensors[f'y{i}'] != expected or mapped != expected:
            raise AssertionError('symbolic polynomial certificate failed')
        certificate.append(dict(output=f'y{i}', surviving_monomials=[f'{a}*{b}' for a,b in sorted(expected)], verified=True))
    matrix = [[column[(i-j)%4] if circular else (column[i-j] if i>=j else row[j-i]) for j in range(4)] for i in range(4)]
    # Independent scalar multiplication: carryless product + polynomial division.
    def multiply(a,b):
        result = reduce(xor, (a<<k for k in range(8) if b>>k&1), 0)
        for k in range(14,7,-1):
            if result>>k&1:
                result ^= 0x11B << (k-8)
        return result
    expected = [reduce(xor,(multiply(a,b) for a,b in zip(r,vector)),0) for r in matrix]
    outputs = [values[f'y{i}'] for i in range(4)]
    if outputs != expected or tuple(outputs) != iterative_result:
        raise AssertionError('numeric comparison failed')
    trace, peak = _lifetime(steps, {f'b{i}' for i in range(4)}, persistent)
    tally = Counter(step['op'] for step in trace)
    return dict(kind=kind, matrix=matrix, column=column, row=row, vector=vector,
                coefficient_inputs=dict(zip(input_names,raw_values)), matrix_steps=coeff_steps,
                persistent_coefficients=[dict(id=name,value=values[name]) for name in sorted(persistent)],
                products=products, formula_terms=term_map, vector_steps=trace, certificate=certificate,
                output=outputs, expected=expected,
                counts=dict(coefficient_xors=len(coeff_steps), vector_xors=tally['xor'], multiplications=tally['multiply'],
                            identity_copies=tally['copy'], retained_coefficients=len(persistent),
                            peak_slots=peak, schedule_steps=len(trace)))
