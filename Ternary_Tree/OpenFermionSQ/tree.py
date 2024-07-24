from __future__ import annotations

"""Tree transform on fermionic operators."""

from openfermion.ops.operators import (FermionOperator, MajoranaOperator,
                                       QubitOperator)
from openfermion.ops.representations import InteractionOperator
from openfermion.utils.operator_utils import count_qubits

from ..Tree.TernaryTree import TernaryTree


def tree(operator: FermionOperator | MajoranaOperator | InteractionOperator,
         n_qubits: int = None,
         tt: TernaryTree = None):
    """Apply the Tree transform.

    Implementation from arXiv:quant-ph/0003137 and
    "A New Data Structure for Cumulative Frequency Tables" by Peter M. Fenwick.

    Note that this implementation is equivalent to the one described in
    arXiv:1208.5986, and is different from the one described in
    arXiv:1701.07072. The one described in arXiv:1701.07072 is implemented
    in OpenFermion as `tree_tree`.

    Args:
        operator (openfermion.ops.FermionOperator):
            A FermionOperator to transform.
        n_qubits (int|None):
            number of qubits
        tt (TernaryTree|None):
            Ternary Tree to build ladder (majorana fermions) operators

    Returns:
        transformed_operator: An instance of the QubitOperator class.

    Raises:
        ValueError: Invalid number of qubits specified.
    """
    if tt is None:
        tt = TernaryTree()
        tt.build_full_tt(n_qubits)
    if isinstance(operator, FermionOperator):
        return _tree_fermion_operator(operator, tt)
    if isinstance(operator, MajoranaOperator):
        return _tree_majorana_operator(operator, tt)
    # if isinstance(operator, InteractionOperator):
    #     return _tree_interaction_operator(operator, tt)
    raise TypeError("Couldn't apply the Tree Transform to object "
                    "of type {}.".format(type(operator)))


def _tree_majorana_operator(operator, tt):
    # Compute the number of qubits.
    N = count_qubits(operator)
    n_qubits = tt.num_nodes

    if n_qubits < N:
        raise ValueError('Invalid number of qubits specified.')

    # Compute transformed operator.
    transformed_terms = (_transform_majorana_term(term=term,
                                                  coefficient=coeff,
                                                  tt=tt)
                         for term, coeff in operator.terms.items())
    return inline_sum(summands=transformed_terms, seed=QubitOperator())


def _transform_majorana_term(term, coefficient, tt):
    # Build the Tree transformed operators.
    transformed_ops = (_transform_majorana_operator(majorana_index, tt)
                       for majorana_index in term)
    return inline_product(factors=transformed_ops,
                          seed=QubitOperator((), coefficient))


def _transform_majorana_operator(majorana_index, tt):
    # q, b = divmod(majorana_index, 2)

    return QubitOperator(tt.get_majorana(majorana_index))


def _transform_operator_term(term: list[tuple[int, int]],
                             coefficient: float,
                             tt: TernaryTree):
    """
    Args:
        term (list[tuple[int, int]]):
            A list of (mode, raising-vs-lowering) ladder operator terms.
        coefficient (float):
        tt (TernaryTree):
    Returns:
        QubitOperator:
    """

    # Build the Tree transformed operators.
    transformed_ladder_ops = (_transform_ladder_operator(
        ladder_operator, tt) for ladder_operator in term)
    return inline_product(factors=transformed_ladder_ops,
                          seed=QubitOperator((), coefficient))


def _tree_fermion_operator(operator, tt):
    # Compute the number of qubits.
    N = count_qubits(operator)
    n_qubits = tt.num_nodes
    if n_qubits < N:
        raise ValueError('Invalid number of qubits specified.')

    # Compute transformed operator.
    transformed_terms = (_transform_operator_term(
        term=term, coefficient=operator.terms[term], tt=tt)
                         for term in operator.terms)
    return inline_sum(summands=transformed_terms, seed=QubitOperator())


#TODO
def _transform_ladder_operator(ladder_operator, tt):
    """
    Args:
        ladder_operator (tuple[int, int]): the ladder operator
        n_qubits (int): the number of qubits
    Returns:
        QubitOperator
    """
    index, action = ladder_operator

    # Initialize the transformed majorana operator (a_p^\dagger + a_p) / 2
    transformed_operator = QubitOperator(tt.get_majorana(2*index + 1), .5)
    # Get the transformed (a_p^\dagger - a_p) / 2
    # Below is equivalent to X(update_set) * Z(parity_set ^ occupation_set)
    transformed_majorana_difference = QubitOperator(
        tt.get_majorana(2 * index + 2), -.5j)

    # Raising
    if action == 1:
        transformed_operator += transformed_majorana_difference
    # Lowering
    else:
        transformed_operator -= transformed_majorana_difference

    return transformed_operator


def inline_sum(summands, seed):
    """Computes a sum, using the __iadd__ operator.
    Args:
        seed (T): The starting total. The zero value.
        summands (iterable[T]): Values to add (with +=) into the total.
    Returns:
        T: The result of adding all the factors into the zero value.
    """
    for r in summands:
        seed += r
    return seed


def inline_product(factors, seed):
    """Computes a product, using the __imul__ operator.
    Args:
        seed (T): The starting total. The unit value.
        factors (iterable[T]): Values to multiply (with *=) into the total.
    Returns:
        T: The result of multiplying all the factors into the unit value.
    """
    for r in factors:
        seed *= r
    return seed


# def _tree_interaction_operator(interaction_operator, n_qubits):
#     """Implementation of the Tree transformation for OpenFermion
#     Interaction Operators. This implementation is equivalent to that described
#     in arXiv:1208.5986, and has been written to optimize compute time by using
#     algebraic expressions for general products a_i^dagger a_j^dagger as outlined
#     in Table II of Seeley, Richard, Love.
#     """
#
#     one_body = interaction_operator.one_body_tensor
#     two_body = interaction_operator.two_body_tensor
#     constant_term = interaction_operator.constant
#
#     # Compute the number of qubits.
#     N = len(one_body)
#     if n_qubits is None:
#         n_qubits = N
#     if n_qubits < N:
#         raise ValueError('Invalid number of qubits specified.')
#
#     qubit_hamiltonian = QubitOperator()
#     qubit_hamiltonian_op = []
#     qubit_hamiltonian_coef = []
#
#     #  For cases A - F see Table II of Seeley, Richard, Love
#     for i in range(n_qubits):
#         # A. Number operators: n_i
#         if abs(one_body[i, i]) > 0:
#             qubit_hamiltonian += _qubit_operator_creation(
#                 *_seeley_richard_love(i, i, one_body[i, i], n_qubits))
#
#         for j in range(i):
#             # Case B: Coulomb and exchange operators
#             if abs(one_body[i, j]) > 0:
#                 operators, coef_list = _seeley_richard_love(
#                     i, j, one_body[i, j], n_qubits)
#                 qubit_hamiltonian_op.extend(operators)
#                 qubit_hamiltonian_coef.extend(coef_list)
#
#                 operators, coef_list = _seeley_richard_love(
#                     j, i, one_body[i, j].conj(), n_qubits)
#                 qubit_hamiltonian_op.extend(operators)
#                 qubit_hamiltonian_coef.extend(coef_list)
#
#             coef = _two_body_coef(two_body, i, j, j, i) / 4
#             if abs(coef) > 0:
#                 qubit_hamiltonian_op.append(
#                     tuple((index, "Z") for index in _occupation_set(i)))
#                 qubit_hamiltonian_op.append(
#                     tuple((index, "Z") for index in _occupation_set(j)))
#                 qubit_hamiltonian_op.append(
#                     tuple((index, "Z") for index in _F_ij_set(i, j)))
#
#                 qubit_hamiltonian_coef.append(-coef)
#                 qubit_hamiltonian_coef.append(-coef)
#                 qubit_hamiltonian_coef.append(coef)
#                 constant_term += coef
#
#     # C. Number-excitation operators: n_i a_j^d a_k
#     for i in range(n_qubits):
#         for j in range(n_qubits):
#             for k in range(j):
#                 if i not in (j, k):
#                     coef = _two_body_coef(two_body, i, j, k, i)
#
#                     if abs(coef) > 0:
#                         number = _qubit_operator_creation(
#                             *_seeley_richard_love(i, i, 1, n_qubits))
#
#                         excitation_op, excitation_coef = _seeley_richard_love(
#                             j, k, coef, n_qubits)
#                         operators_hc, coef_list_hc = _seeley_richard_love(
#                             k, j, coef.conj(), n_qubits)
#                         excitation_op.extend(operators_hc)
#                         excitation_coef.extend(coef_list_hc)
#                         excitation = _qubit_operator_creation(
#                             excitation_op, excitation_coef)
#
#                         number *= excitation
#                         qubit_hamiltonian += number
#
#     # D. Double-excitation operators: c_i^d c_j^d c_k c_l
#     for i in range(n_qubits):
#         for j in range(i):
#             for k in range(j):
#                 for l in range(k):
#                     coef = -_two_body_coef(two_body, i, j, k, l)
#                     if abs(coef) > 0:
#                         qubit_hamiltonian += _hermitian_one_body_product(
#                             i, j, k, l, coef, n_qubits)
#
#                     coef = -_two_body_coef(two_body, i, k, j, l)
#                     if abs(coef) > 0:
#                         qubit_hamiltonian += _hermitian_one_body_product(
#                             i, k, j, l, coef, n_qubits)
#
#                     coef = -_two_body_coef(two_body, i, l, j, k)
#                     if abs(coef) > 0:
#                         qubit_hamiltonian += _hermitian_one_body_product(
#                             i, l, j, k, coef, n_qubits)
#
#     qubit_hamiltonian_op.append(())
#     qubit_hamiltonian_coef.append(constant_term)
#     qubit_hamiltonian += _qubit_operator_creation(qubit_hamiltonian_op,
#                                                   qubit_hamiltonian_coef)
#
#     return qubit_hamiltonian
#
#
# def _two_body_coef(two_body, a, b, c, d):
#     return two_body[a, b, c, d] - two_body[a, b, d, c] + two_body[
#         b, a, d, c] - two_body[b, a, c, d]
#
#
# def _hermitian_one_body_product(a, b, c, d, coef, n_qubits):
#     """ Takes the 4 indices for a two-body operator and constructs the
#     Tree form by splitting the two-body into 2 one-body operators,
#     multiplying them together and then re-adding the Hermitian conjugate to
#     give a Hermitian operator. """
#
#     c_dag_c_ac = _qubit_operator_creation(
#         *_seeley_richard_love(a, c, coef, n_qubits))
#     c_dag_c_bd = _qubit_operator_creation(
#         *_seeley_richard_love(b, d, 1, n_qubits))
#     c_dag_c_ac *= c_dag_c_bd
#     hermitian_sum = c_dag_c_ac
#
#     c_dag_c_ca = _qubit_operator_creation(
#         *_seeley_richard_love(c, a, coef.conj(), n_qubits))
#     c_dag_c_db = _qubit_operator_creation(
#         *_seeley_richard_love(d, b, 1, n_qubits))
#     c_dag_c_ca *= c_dag_c_db
#
#     hermitian_sum += c_dag_c_ca
#     return hermitian_sum
#
#
# def _qubit_operator_creation(operators, coefficents):
#     """ Takes a list of tuples for operators/indices, and another for
#     coefficents"""
#
#     qubit_operator = QubitOperator()
#
#     for index in zip(operators, coefficents):
#         qubit_operator += QubitOperator(index[0], index[1])
#
#     return qubit_operator
#
#


