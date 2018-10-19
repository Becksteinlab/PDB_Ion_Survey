import pdbionsurvey.coordination as coordination

def test_get_charge():
    assert(all([coordination.get_charge(i) == 1 for i in ['LI', 'NA', 'K', 'RB', 'CS', 'TL', 'RH', 'AG', 'AU']]))
    assert(all([coordination.get_charge(i) == 2 for i in ['MG', 'CA', 'SR', 'BA', 'MN', 'CO', 'NI', 'PD', 'PT', 'CU', 'ZN', 'CD', 'HG', 'PB']]))
    assert(all([coordination.get_charge(i) == 3 for i in ['LA',  'V', 'CR', 'FE', 'RU', 'OS', 'AL', 'GA', 'IN', 'SB']]))
    assert(all([coordination.get_charge(i) == 4 for i in ['ZR', 'IR']]))
    assert(all([coordination.get_charge(i) == 6 for i in ['W']]))
    assert(all([coordination.get_charge(i) == -1 for i in ['F', 'CL', 'BR', 'IOD', 'I']]))

