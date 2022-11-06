#!/home/mario/Projects/boolean_2/software/venv/bin/venv python

# Libraries
import itertools
from random import choice, sample
from string import ascii_uppercase, digits
from sympy import SOPform
from quine_mccluskey.qm import QuineMcCluskey
# from pyboolnet.Attractors import compute_attractors_tarjan
# from pyboolnet.FileExchange import bnet2primes
# from pyboolnet.StateTransitionGraphs import primes2stg
# from exceptions import NoSolutionException, InputModificationException


# Functions
def left_zfill(word, n_digits):
    """
    DESCRIPTION:
    Given a string a number, the function introduces zeros by the left side 
    until the string the size specified in the number.
    :param word: [string] the string that is to be extended. 
    :param n_digits: [int] the number of digits of the final word. 
    :return: [string] extended word. 
    """
    return '0' * (n_digits - len(word)) + word
    

# def function_simplifier(node, nodes, n_nodes, terms, method='Quine-McCluskey'):
#     """
#     DESCRIPTION:
#     A function to simplify the terms of a given boolean function.
#     :param node: [str] node that denotes the boolean function that we
#     are simplifying.
#     :param nodes: [list] variables of the function to simplify.
#     :param n_nodes: [int] the number of nodes computed to save time.
#     :param terms: [frozenset] strings the represent the terms.
#     :param method: [str] variable to indicate the simplification method to
#     use.
#     :return: [list] minterms of the simplfied function.
#     """
#     if not terms:
#         raise NoSolutionException()
#     if method == 'Quine-McCluskey':
#         qm = QuineMcCluskey(use_xor=False)
#         minterms = [int(term, 2) for term in terms]
#         # Special cases
#         if minterms == [0]:
#             results = [left_zfill(str(minterm), n_nodes) for minterm in minterms]
#         # General case
#         else:
#             results = [left_zfill(minterm.replace('-', '*'), n_nodes) for minterm
#                 in qm.simplify(minterms)]
#     elif method == 'SymPy':
#         terms = [[int(digit) for digit in list(term)] for term in terms]
#         simplified = str(SOPform(nodes, terms))
#         terms = [
#             term.replace('(', '').replace(')', '').replace(' ', '').split('&')
#             for term in simplified.split('|')
#         ]
#         results = []
#         for term in terms:
#             result = dict(zip(nodes, ['*'] * n_nodes))
#             for factor in term:
#                 result[factor.replace('~', '')] = '0' if '~' in factor else '1'
#             results.append(''.join(result.values()))  
#     else:
#         raise AttributeError('Introduced non-valid simplification mode')
#     return results

def minterms2bnet(variables, minterms):
    """
    DESCRIPTION:
    A function to obtain the function expression from the minterms in boolnet format.
    It is supposed that all the nodes are in alphabetical order. Important, this 
    function only returns the SOP without the name of the function, what is at the 
    left side of the comma.
    :param nodes: [tuple] the function variables.
    :param minterms: [frozenset] the minterms to build the expression.
    :return: [str] the function expression in boolnet format.
    """
    # When the function is always false
    if not minterms:
        return '0'
    # Check when the minterms are the whole space
    if len(minterms) == 2 ** len(variables):
        return '1'
    # Generañ case
    qm = QuineMcCluskey(use_xor=False)
    simplified_minterms = qm.simplify([int(term, 2) for term in minterms])
    simplified_expression = [left_zfill(minterm.replace('-', '*'), len(variables)) for minterm in simplified_minterms]
    # Pass the expression to boolnet format
    n_variables = range(len(variables))
    bnet_expression = ' | '.join([
        '&'.join([variables[i] if int(minterm[i]) else '!' + variables[i] for i in n_variables if minterm[i] != '*'])
        for minterm in simplified_expression
    ])
    return bnet_expression


# def network_formatter(network, min_attractors=2, max_attractors=4):
#     """
#     DESCRIPTION:
#     A function to compute the attractors and format a given network to prepare it for 
#     the storage.
#     :param network: [dict] the boolean network whose attractors are going to be 
#     computed.
#     :param min/max_attractors: [int] a filtering parameters. The only networks of interest
#     are those with a number of attractors within the limits.
#     :return: [dict] the boolean network with the attractors stored. 
#     """
#     # Write every net function in BoolNet format
#     network['network_terms'] = network['network'].copy()
#     network['network'] = {
#         node: f"{node}, {minterms2bnet(network['nodes'], network['network'][node])}" 
#         for node in network['nodes']
#     }
#     # Obtain the state transition graph
#     primes = bnet2primes('\n'.join(network['network'].values()))
#     stg = primes2stg(primes, "synchronous")
#     # Obtain the attractors
#     steady, cyclic = compute_attractors_tarjan(stg)
#     network['attractors'] = {'steady': steady, 'cyclic': cyclic}
#     # Return only if the attractor condition is met
#     if min_attractors <= len(steady + cyclic) and len(steady + cyclic) <= max_attractors:
#         # Obtain the initial expression
#         network['pre_network'] = {
#             node: f"{node}, {minterms2bnet(network['nodes'], network['pre_network'][node])}" 
#             for node in network['nodes']
#         }
#         # Simplify both initial and final expressions
#         # Delete not needed variables
#         del network['max_iterations']
#         del network['graph_space']
#         del network['nodes']
#         del network['priority']
#         del network['network_terms']
#         # Return
#         return network


# def filter_boolean_networks(boolean_networks, attractors=None, n_attractors=None, partial=False):
#     """
#     DESCRIPTION:
#     A function to filter the boolean networks according to different 
#     criteria based on the attractors.
#     :param boolean_networks: [list] the boolean networks to filter.
#     :param attractors: [dict] the attractors that we want the networks
#     to show. Format: {'steady': [...], 'cyclic': [...]}
#     :param n_attractors: [int] the number of attractors we want the
#     networks to show, steady + cyclic.
#     :param partial: [bool] if partial == True the network showing at
#     least one of the specified attractors passes the filter, also the
#     network having a number of attractors <= n_attractors. If 
#     partial == False, the network should show all the attractors and
#     have a number of attractors equal to n_attractors.
#     :return: [list] the boolean networks that meet the criteria.
#     """
#     # Kernels
#     def by_attractor(network):
#         """
#         DESCRIPTION:
#         A function to filter by attractor. It is returned true only if
#         the network has the attractors specified in attractors.
#         :param network: [dict] the boolean network to filter.
#         :return: [bool] value to indicate if the network meet the 
#         criterium.
#         """
#         # Obtain all the attractors from the network
#         network_attractors = network['attractors']['steady'] +\
#             network['attractors']['cyclic']
#         # If partial, look only for one attractor.
#         if partial:
#             condition = any([attractor in network_attractors 
#                 for attractor in attractors])
#         else:
#             condition = all([attractor in network_attractors 
#                 for attractor in attractors])
#         return condition

#     def by_n_attractors(network):
#         """
#         DESCRIPTION:
#         A function to filter by the number of attractors specified.
#         :param network: [dict] the boolean network to filter.
#         :return: [bool] value to indicate if the network meet the 
#         criterium.
#         """
#         # Obtain the number of attractors in the network
#         n_network_attractors = len(network['attractors']['steady'] +\
#             network['attractors']['cyclic'])
#         # If partial, look only for networks that at most have the 
#         # specified number of attractors
#         condition = (n_network_attractors <= n_attractors 
#             if partial else n_network_attractors == n_attractors)
#         return condition

#     # Filtering
#     if attractors is not None:
#         boolean_networks = list(filter(by_attractor, boolean_networks))
#     if n_attractors is not None:
#         boolean_networks = list(filter(by_n_attractors, boolean_networks))
#     return boolean_networks


def prefilter_by_attractor(networks, attractors):
    """
    DESCRIPTION:
    A generator to filter boolean networks based on if they hold certain
    attractors or not.
    :param networks: [list] boolean networks (dict) to filter based on
    their attractors.
    :param attractors: [list] attractors (str) to filter the networks.
    :return: [dict] network that shows all the attractors.
    """
    # Helper functions
    def check_node(node_index, nodes, attractor, network):
        node = nodes[node_index]
        # Check if the attractor value is 1 or 0
        if int(attractor[node_index]):
            result = attractor in network[node]
        else:
            result = attractor not in network[node]
        return result

    # Iterate over every network
    nodes = list(networks[0].keys())
    for network in networks:
        # Network might be a None resulting from unsuccessful conflict solving
        if network is not None:
            attractor_conditions = []
            # Check every attractor
            for attractor in attractors:
                node_conditions = [check_node(node_index, nodes, attractor, network) 
                    for node_index in range(len(nodes))]
                attractor_conditions.append(all(node_conditions))
            if all(attractor_conditions):
                yield network