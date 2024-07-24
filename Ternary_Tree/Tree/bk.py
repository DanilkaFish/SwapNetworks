def _update_set(index, n_qubits):
    """The bits that need to be updated upon flipping the occupancy
    of a mode."""
    indices = set()

    # For bit manipulation we need to count from 1 rather than 0
    index += 1
    # Ensure index is not a member of the set
    index += index & -index

    while index <= n_qubits:
        indices.add(index - 1)
        # Add least significant one to index
        # E.g. 00010100 -> 00011000
        index += index & -index
    return indices


def _occupation_set(index):
    """The bits whose parity stores the occupation of mode `index`."""
    indices = set()

    # For bit manipulation we need to count from 1 rather than 0
    index += 1

    indices.add(index - 1)
    parent = index & (index - 1)
    index -= 1
    while index != parent:
        indices.add(index - 1)
        # Remove least significant one from index
        # E.g. 00010100 -> 00010000
        index &= index - 1
    return indices


def _parity_set(index):
    """The bits whose parity stores the parity of the bits 0 .. `index`."""
    indices = set()

    while index > 0:
        indices.add(index - 1)
        # Remove least significant one from index
        # E.g. 00010100 -> 00010000
        index &= index - 1
    return indices


def transform_majorana_operator(majorana_index, n_qubits):
    q, b = divmod(majorana_index, 2)

    update_set = _update_set(q, n_qubits)
    update_set.add(q)
    occupation_set = _occupation_set(q)
    parity_set = _parity_set(q)

    if b:
        return ([(q, 'Y')] + [(i, 'X') for i in update_set - {q}] +
                             [(i, 'Z')
                              for i in (parity_set ^ occupation_set) - {q}])
    else:
        return ([(i, 'X') for i in update_set] +
                             [(i, 'Z') for i in parity_set])


def bk_majorana_operators(n_qubits):
    branches = []
    for i in range(2*n_qubits):
        branches = branches + [transform_majorana_operator(i, n_qubits)]
    return branches
