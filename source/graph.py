
# Libraries
import os
import itertools
from random import sample
from string import ascii_uppercase, digits
from source.ncbf_utils import ncbf_generator
from source.bn_utils import prefilter_by_attractor
from source.bn_utils import minterms2bnet


# Classes
class Graph:

    # Methods
    def __init__(self, activators, inhibitors, attractors, networks_path):
        """
        DESCRIPTION:
        The constructor of the Graph object. All the network inference
        is based on this object.
        :param activators: [dict] dict in which the keys are the nodes, and the
        values their activators.
        :param inhibitors: [dict] the same but for the inhibitors.
        :param attractors: [list] a list with the searched attractors in str 
        format.
        :param networks_path: [str] path to folder to print the networks in.
        """
        # Always, the nodes are ordered alphabetically
        self.nodes = tuple(sorted(activators.keys()))
        self.activators = {key: tuple(sorted(activators[key])) 
            for key in activators.keys()}
        self.inhibitors = {key: tuple(sorted(inhibitors[key])) 
            for key in inhibitors.keys()}
        self.n_nodes = len(self.nodes)
        self.attractors = None if attractors == "None" else attractors
        # Generate all the possible minterms in a space of len(nodes) variables.
        # IMPORTANT: the node position in every term is alphabetical: A:0, B:1...
        self.graph_space = frozenset(
            '{:0{}b}'.format(i, self.n_nodes) for i in range(2 ** self.n_nodes)
            )
        # Check for input nodes
        self.input_nodes = tuple([node for node in self.nodes 
            if not self.activators[node] and not self.inhibitors[node]])
        # Path to print the networks
        self.networks_path = networks_path

    def __str__(self):
        """
        DESCRIPTION:
        Method to create the object string representation.
        :return: [str] object string representation.
        """
        representation = f"""
        ************************************************************************
        Graph object
        ************************************************************************
        Nodes:
        {self.nodes}
        Activators:
        {self.activators}
        Inhibitors:
        {self.inhibitors}
        N_nodes:
        {self.n_nodes}
        Attractors:
        {self.attractors}
        Input_nodes:
        {self.input_nodes}
        Number of computed pathway groups:
        {len(self.get_pathway_groups())}
        Number of computed nested canalised boolean functions:
        {len(self.get_ncbf_networks())}
        Number of computed nested canalised boolean functions after filtering:
        {len(self.get_filtered_ncbf_networks())}
        Graph space size:
        {len(self.graph_space)}
        ************************************************************************
        """
        return representation

    def get_pathway_groups(self):
        """
        DESCRIPTION:
        """
        if "pathway_groups" in dir(self):
            return self.pathway_groups
        else:
            return []

    def get_ncbf_networks(self):
        """
        DESCRIPTION:
        """
        if "ncbf_networks" in dir(self):
            return self.ncbf_networks
        else:
            return []
    
    def get_filtered_ncbf_networks(self):
        """
        DESCRIPTION:
        """
        if "ncbf_networks" in dir(self):
            return self.filtered_ncbf_networks
        else:
            return []
        
    def get_networks(self):
        """
        DESCRIPTION:
        """
        if "networks" in dir(self):
            return self.networks
        else:
            return []

    def obtain_pathways_from_graph(self):
        """
        DESCRIPTION:
        The method to obtain all the possible groups of pathways from the graph 
        based on the pairs of canalising/canalised values described in the 
        article. The groups are added to the Graph object. Every pathway is a 
        tuple in which the fields have the following meaning:
        1. Antecedent: [str] expression describing the cause of the effect.
        2. Consequent: [str] expression describing the nodes that suffers the 
        effect.
        3. Activator: [bool] a flag to indicate if the pathway is an activator 
        of the consequent (True) or not (False).
        4. Domain: [set] strings that represent the minterms of the expression
        present in the left side of the pathway. The right side is always
        exactly the letter shown in the consequent field, there is no need to 
        represent that function.
        """
        # Helper functions
        def pathway_serializer(antecedent, consequent, definition):
            """
            DESCRIPTION: 
            Helper function for the list comprehension of below. 
            :param antecedent: [str] value to write in the left side of the
            pathway.
            :param consequent: [str] value to write in the right side of the
            pathway.
            :param definition: [tuple] canalising/canalised value pair of the
            pathway.
            :return: [dict] formatted pathway.
            """
            return {
                'id': ''.join(sample(ascii_uppercase + digits, 5)),
                'antecedent': antecedent,
                'consequent': consequent,
                'activator': bool(definition[1]),
                'domain': frozenset(filter(
                    lambda term: term[self.nodes.index(antecedent)] == str(definition[0]),
                    self.graph_space))}
        
        def pathway_manager(pathway_group, input_pathways):
            """
            DESCRIPTION:
            Helper function to organise the pathways first by node and after
            by activators and inhibitors. In the code, every activator will have
            a canalised value of 1 and vice versa. The function changes the sense
            of activator and inhibitor. From now on, activator has a canalised
            value of 1 and inhibitor of 0.
            :param pathway_group: [tuple] the two lists of activators and 
            inhibitors.
            :param input_pathways: [dict] pathways of the input nodes. They are 
            the same for every group.
            :return: [dict] the pathways grouped by node and activators/inhibitors.
            """
            # Format the pathways
            pathways = [it for sb in pathway_group for it in sb]
            pathways = {node: {'activators': list(filter(lambda pathway: (pathway['consequent'] == node) and pathway['activator'], pathways)),
                                'inhibitors': list(filter(lambda pathway: (pathway['consequent'] == node) and not pathway['activator'], pathways))}
                for node in self.nodes}
            # Introduce the input pathways
            [pathways.update({node: input_pathways[node]}) for node in input_pathways.keys()]
            return pathways
        
        # Create all the pathways with both canalising/canalised pairs
        activator_pathways = [[None] * len(self.activators[node]) for node in self.nodes]
        inhibitor_pathways = [[None] * len(self.inhibitors[node]) for node in self.nodes]
        i = 0
        for node in self.nodes:
            j = 0
            for activator in self.activators[node]:
                activator_pathways[i][j] = [
                    pathway_serializer(activator, node, (0, 0)),
                    pathway_serializer(activator, node, (1, 1))
                ]
                j += 1
            j = 0
            for inhibitor in self.inhibitors[node]:
                inhibitor_pathways[i][j] = [
                    pathway_serializer(inhibitor, node, (0, 1)),
                    pathway_serializer(inhibitor, node, (1, 0))
                ]
                j += 1
            i += 1
        activator_pathways = [it for sb in activator_pathways for it in sb]
        activator_pathways = itertools.product(*activator_pathways)
        inhibitor_pathways = [it for sb in inhibitor_pathways for it in sb]
        inhibitor_pathways = itertools.product(*inhibitor_pathways)
        pathways = itertools.product(activator_pathways, inhibitor_pathways)
        # Obtain the pathways of the inputs (from and to themselves)
        input_pathways = {node: [] for node in self.input_nodes}
        for node in self.input_nodes:
            input_pathways[node] = {
                'activators': [pathway_serializer(node, node, (1, 1))],
                'inhibitors': [pathway_serializer(node, node, (0, 0))]
            }
        # Organise every group by node first and by activator/inhibitor second
        self.pathway_groups = [pathway_manager(group, dict(input_pathways)) for group in pathways]

    def generate_priority_matrices(self):
        """
        DESCRIPTION:
        A method that adds to the Graph object the priority matrices used in the 
        inference. There are two priority matrices per simulation, one for activators
        and another for inhibitors.
        """
        # Obtain all the possible node combinations
        combinations = [itertools.combinations(self.nodes, i + 1)
            for i in range(self.n_nodes)]
        combinations = sorted([''.join(it) for sb in combinations for it in sb])
        # Generate all the possible priorities
        # We set arbitrarily all the priorities between 0 and 1000, where 0 is the
        # lowest priority and 1000 the highest
        priorities = range(0, 1000)
        # Activators
        activator_matrices = []
        # We generate as many priority matrices as simulations that we will perform
        for _ in range(self.n_simulations):
            # Generate the random matrices
            activator_matrices.append({
                key: dict(zip(combinations, sample(priorities, len(combinations)))) 
                for key in combinations
                })
        # Inhibitors
        inhibitor_matrices = []
        for _ in range(self.n_simulations):
            # Generate the random matrices
            inhibitor_matrices.append({
                key: dict(zip(combinations, sample(priorities, len(combinations)))) 
                for key in combinations
                })
        # Store the matrices
        self.priority_matrices = zip(activator_matrices, inhibitor_matrices)
    
    def generate_NCBFs(self):
        """
        DESCRIPTION:
        A method to build all the groups of NCBF based on the pathway groups. 
        Multiple NCBF groups are built per pathway group. The information about the 
        canalising/canalised pairs is already stored in the pathways. The NCBF
        groups are added to the Graph object.
        """
        # Helper functions
        # There two depending on whether the pathways groups are mixed or not for
        # the ncbf and the conflicts strategy
        def ncbf_formatter_standard(ncbf_group, pathway_group_id):
            # This function assigns the ID of the pathways used to build the network
            networks = [dict(zip(self.nodes, network)) for network in itertools.product(*ncbf_group)]
            return zip([pathway_group_id] * len(networks), networks)

        def ncbf_formatter_mixed_pathways(ncbf_group):
            # This assigns the prefixed ID of pathways independently of the network
            networks = [dict(zip(self.nodes, network)) for network in itertools.product(*ncbf_group)]
            return zip([self.pathways_index] * len(networks), networks)

        # Generate by node all the NCBFs
        total_ncbf = [[ncbf_generator(group[node]['activators'], group[node]['inhibitors'], self.graph_space, set(self.nodes)) 
            for node in self.nodes] for group in self.pathway_groups]
        # # Format all the NCBF groups conveniently: (pathway group position, NCBF
        # # network in dict). 
        # if self.mixed_pathways:
        #     pre_networks = [it for sb in [ncbf_formatter_mixed_pathways(total_ncbf[i]) for i in range(len(total_ncbf))] for it in sb]
        # else:
        #     pre_networks = [it for sb in [ncbf_formatter_standard(total_ncbf[i], i) for i in range(len(total_ncbf))] for it in sb]
        ncbf_networks = [it for sb in [ncbf_formatter_standard(total_ncbf[i], i) for i in range(len(total_ncbf))] for it in sb]
        # Filter the equivalent networks
        codes = []
        final_ncbf_networks = []
        for network in ncbf_networks:
            # code = str(network[0]) + '$$' + '&'.join(['|'.join(sorted(net)) for _, net in sorted(network[1].items(), key=lambda x: x[0])])
            code = '&'.join(['|'.join(sorted(net)) for _, net in sorted(network[1].items(), key=lambda x: x[0])])
            if code not in codes:
                final_ncbf_networks.append(network)
                codes.append(code)
        # Store NCBF networks
        # NOTE: remove ID because we do not need to match with the pathways again.
        self.ncbf_networks = [net[1] for net in final_ncbf_networks]
        self.networks = self.ncbf_networks
    
    def prefilter(self):
        """
        DESCRIPTION:
        A method to select only the networks among which the searched attractors are present.
        A cheap prefiltering before computing the attractors with the Tarjan algorithm and 
        PyBoolNet.
        """
        if (self.attractors is not None) and (self.attractors != []):
            self.filtered_ncbf_networks = list(filter(lambda network: network is not None, self.get_ncbf_networks()))
            if self.filtered_ncbf_networks:
                self.filtered_ncbf_networks = list(prefilter_by_attractor(self.filtered_ncbf_networks, self.attractors))
            self.networks = self.filtered_ncbf_networks
            print(f'Total networks after prefiltering: {len(self.filtered_ncbf_networks)}')
        else:
            print('No attractor-based prefiltering performed')
            self.filtered_ncbf_networks = self.get_ncbf_networks()
            self.networks = self.filtered_ncbf_networks

    def print_networks_to_folder(self, folder_path=None, prefix="hallo"):
        """
        DESCRIPTION:
        Method to write the networks into a specified folder.
        :param folder_path: [str] path to the folder to store the networks.
        :param prefix: [str] prefix to name the network files.
        """
        if not folder_path:
            folder_path = self.networks_path
        networks = [
            "\n".join([f"{node}, " + minterms2bnet(self.nodes, network[node]) for node in self.nodes])
            for network in self.get_networks()
        ]
        i = 0
        os.system(f"mkdir -p {folder_path}")
        for network in networks:
            with open(folder_path + f"/{prefix}_{i}.txt", "w") as file:
                file.write(network)
            i += 1
        return None
