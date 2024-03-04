"""This module parses biomedical knowledge bases into several types of 
objects (dictionaries and Networkx graph)."""

import csv
import networkx as nx
import obonet

data_dir = "data/kbs/"


class KnowledgeBase:
    """Represents a knowledge base that is loaded from a given local file."""

    def __init__(
        self,
        kb,
        terms_filename=None,
        edges_filename=None,
        kb_filename=None,
        file_format=None,
    ):
        self.kb = kb
        self.terms_filename = terms_filename
        self.edges_filename = edges_filename
        self.kb_filename = kb_filename
        self.file_format = file_format  # obo, tsv, csv, txt

        self.root_dict = {
            "go_bp": ("GO:0008150", "biological_process"),
            "go_cc": ("GO:0005575", "cellular_component"),
            "chebi": ("CHEBI:00", "root"),
            "hp": ("HP:0000001", "All"),
            "medic": ("MESH:C", "Diseases"),
            "ctd_anat": ("MESH:A", "Anatomy"),
            "ctd_chem": ("MESH:D", "Chemicals"),
            "ctd_gene": ("", ""),
            "do": "DOID:4",
        }

        self.name_to_id = None
        self.id_to_info = None
        self.synonym_to_id = None
        self.edges = None
        self.graph = None
        self.child_to_parent = None
        self.alt_id_to_id = None
        self.umls_to_hp = None

        # ---------------------------------------------------------------------
        #                 Load the info about the given KB
        # ---------------------------------------------------------------------
        if file_format == "obo":
            self.load_obo()

        elif self.file_format == "tsv":
            self.load_tsv()

        elif self.kb == "ncbi_taxon":
            self.load_ncbi_taxon()

        elif self.kb == "ncbi_gene":
            self.load_ncbi_gene()

        elif self.file_format == "txt":
            self.load_txt()

    def load_obo(self):
        """Load KBs from local .obo files (ChEBI, HP, MEDIC, GO) into
        structured dicts containing the mappings name_to_id, id_to_name,
        id_to_info, synonym_to_id, child_to_parent, umls_to_hp and the list of
        edges between concepts. For 'chebi', only the concepts in the subset
        3_STAR are included, which correpond to manually validated entries.
        """

        filepath = data_dir
        filepaths = {
            "medic": "CTD_diseases",
            "chebi": "chebi",
            "go_bp": "go-basic",
            "go_cc": "go-basic",
            "do": "doid",
            "hp": "hp",
            "cellosaurus": "cellosaurus",
            "cl": "cl-basic",
            "uberon": "uberon-basic",
        }

        if self.kb in filepaths:
            filepath += filepaths[self.kb] + ".obo"

        else:
            filepath += self.kb_filename

        name_to_id = {}
        id_to_name = {}
        id_to_info = {}
        synonym_to_id = {}
        child_to_parent = {}
        alt_id_to_id = {}
        umls_to_hp = {}

        graph = obonet.read_obo(filepath)
        edges = []

        for node in graph.nodes(data=True):
            add_node = False

            if "name" in node[1]:
                node_id, node_name = node[0], node[1]["name"]

                # node_id = node_id.replace(':', '_')

                if self.kb == "go_bp":
                    # For go_bp, ensure that only Biological Process
                    # concepts are considered

                    if node[1]["namespace"] == "biological_process":
                        name_to_id[node_name] = node_id
                        id_to_name[node_id] = node_name
                        # id_to_info[node_id] = node_name
                        add_node = True

                elif self.kb == "go_cc":
                    if node[1]["namespace"] == "cellular_component":
                        name_to_id[node_name] = node_id
                        id_to_name[node_id] = node_name
                        # id_to_info[node_id] = node_name
                        add_node = True

                elif self.kb == "medic":
                    # if node_id[0:4] != "OMIM":
                    #    # Exclude OMIM concepts #TODO: revise later
                    name_to_id[node_name] = node_id
                    id_to_name[node_id] = node_name
                    add_node = True

                else:
                    name_to_id[node_name] = node_id
                    id_to_name[node_id] = node_name
                    add_node = True

                if "alt_id" in node[1].keys():
                    for alt_id in node[1]["alt_id"]:
                        # alt_id_to_id[alt_id.replace(':', '_')] = node_id
                        alt_id_to_id[alt_id] = node_id

                if "is_obsolete" in node[1].keys() and node[1]["is_obsolete"] is True:
                    add_node = False
                    del name_to_id[node_name]
                    del id_to_name[node_id]

                # Check parents for this node
                if "is_a" in node[1].keys() and add_node:
                    # The root node of the ontology does not
                    # have is_a relationships

                    if len(node[1]["is_a"]) == 1:
                        # Only consider concepts with 1 direct ancestor
                        child_to_parent[node_id] = node[1]["is_a"][0]

                    for parent in node[1]["is_a"]:
                        # To build the edges list, consider all
                        # concepts with at least one ancestor
                        edges.append((node_id, parent))

                if self.kb == "cellosaurus":
                    if "relationship" in node[1].keys() and add_node:
                        relations = node[1]["relationship"]

                        for relation in relations:
                            if relation[:13] == "derived_from ":
                                parent = relation.split("derived_from")[1][1:]
                                edges.append((parent, node_id))

                if "synonym" in node[1].keys() and add_node:
                    # Check for synonyms for node (if they exist)

                    for synonym in node[1]["synonym"]:
                        synonym_name = synonym.split('"')[1]
                        synonym_to_id[synonym_name] = node_id

                if "xref" in node[1].keys() and add_node:
                    if self.kb == "hp":
                        # Map UMLS concepts to HPO concepts
                        for xref in node[1]["xref"]:
                            if xref[:4] == "UMLS":
                                umls_id = xref.strip("UMLS:")
                                umls_to_hp[umls_id] = node_id

        if self.kb in self.root_dict:
            root_concept_name = self.root_dict[self.kb][1]
            root_id = ""

            if root_concept_name not in name_to_id:
                root_id = self.root_dict[self.kb][0]
                name_to_id[root_concept_name] = root_id
                id_to_name[root_id] = root_concept_name

        # ----------------------------------------------------------------------
        # Add misssing edges between the ontology root and
        # sub-ontology root concepts
        if self.kb == "chebi":
            chemical_entity = "CHEBI_24431"
            edges.append((chemical_entity, root_id))
            role = "CHEBI_50906"
            edges.append((role, root_id))
            subatomic_particle = "CHEBI_36342"
            edges.append((subatomic_particle, root_id))
            application = "CHEBI_33232"
            edges.append((application, root_id))

        kb_graph = nx.DiGraph([edge for edge in edges])

        # Build id_to_info (KB-ID: (outdegree, indegree, num_descendants))
        for node in kb_graph.nodes:
            num_descendants = len(nx.descendants(kb_graph, node))

            id_to_info[node] = (
                kb_graph.out_degree(node),
                kb_graph.in_degree(node),
                num_descendants,
            )

        node_to_node = {}

        for edge in edges:
            if edge[0] in node_to_node:
                node_to_node[edge[0]].append(edge[1])

            else:
                node_to_node[edge[0]] = [edge[1]]

            if edge[1] in node_to_node:
                node_to_node[edge[1]].append(edge[0])

            else:
                node_to_node[edge[1]] = [edge[0]]

        self.name_to_id = name_to_id
        self.id_to_name = id_to_name
        self.id_to_info = id_to_info
        self.synonym_to_id = synonym_to_id
        self.edges = edges
        self.graph = kb_graph
        self.child_to_parent = child_to_parent
        self.alt_id_to_id = alt_id_to_id
        self.umls_to_hp = umls_to_hp
        self.node_to_node = node_to_node

    def load_tsv(self):
        """Load KBs from local .tsv files (CTD-Chemicals, CTD-Anatomy)
           into structured dicts containing the mappings name_to_id,
           id_to_info, synonym_to_id, child_to_parent, and the list of edges
           between concepts.
        """

        kb_dict = {
            "ctd_chem": "CTD_chemicals",
            "ctd_anat": "CTD_anatomy",
            "ctd_gene": "CTD_genes",
            "medic": "CTD_diseases",
        }
        filepath = f"{data_dir}{self.kb}/{kb_dict[self.kb]}.tsv"

        name_to_id = {}
        id_to_name = {}
        id_to_info = {}
        synonym_to_id = {}
        child_to_parent = {}
        edges = []

        with open(filepath) as kb_file:
            reader = csv.reader(kb_file, delimiter="\t")
            row_count = int()

            for row in reader:
                row_count += 1

                if row_count >= 30:
                    node_name = row[0]
                    node_id = row[1]  # .replace(':', '_')
                    node_parents = row[4].split("|")
                    synonyms = row[7].split("|")
                    name_to_id[node_name] = node_id
                    id_to_name[node_id] = node_name

                    if len(node_parents) == 1:  #
                        # Only consider concepts with 1 direct ancestor
                        child_to_parent[node_id] = node_parents[0]

                    for synonym in synonyms:
                        synonym_to_id[synonym] = node_id

                    for parent in node_parents:
                        # To build the edges list, consider
                        # all concepts with at least one ancestor

                        edges.append((node_id, parent))

        root_concept_name = self.root_dict[self.kb][1]
        root_concept_id = self.root_dict[self.kb][0]
        name_to_id[root_concept_name] = root_concept_id
        id_to_name[root_concept_id] = root_concept_name

        kb_graph = nx.DiGraph([edge for edge in edges])

        # Build id_to_info (KB-ID: (outdegree, indegree, num_descendants))
        for node in kb_graph.nodes:
            num_descendants = len(nx.descendants(kb_graph, node))

            id_to_info[node] = (
                kb_graph.out_degree(node),
                kb_graph.in_degree(node),
                num_descendants,
            )

        node_to_node = {}

        for edge in edges:
            if edge[0] in node_to_node:
                node_to_node[edge[0]].append(edge[1])

            else:
                node_to_node[edge[0]] = [edge[1]]

            if edge[1] in node_to_node:
                node_to_node[edge[1]].append(edge[0])

            else:
                node_to_node[edge[1]] = [edge[0]]

        self.name_to_id = name_to_id
        self.id_to_name = id_to_name
        self.id_to_info = id_to_info
        self.synonym_to_id = synonym_to_id
        self.edges = edges
        self.graph = kb_graph
        self.child_to_parent = child_to_parent
        self.node_to_node = node_to_node

    def load_ncbi_taxon(self):
        """Load KBs from local .csv files (NCBITaxon) into structured dicts
            containing the mappings name_to_id, id_to_info, synonym_to_id.
        """

        filepath = data_dir

        if self.kb == "ncbi_taxon":
            filepath += "NCBITAXON.csv"

        name_to_id = {}
        id_to_name = {}
        id_to_info = {}
        synonym_to_id = {}
        child_to_parent = {}
        edges = []

        with open(filepath) as csv_file:
            reader = csv.reader(csv_file, delimiter=",")
            row_count = int()

            for row in reader:
                row_count += 1

                if row_count > 1 and "NCBITAXON/" in row[0]:
                    rank_node = row[9]

                    if rank_node == "species":
                        node_name = row[1]
                        node_id = "NCBITaxon_" + row[0].split("NCBITAXON/")[1]
                        synonyms = row[2].split("|")
                        name_to_id[node_name] = node_id
                        id_to_name[node_id] = node_name

                        if row[7] != "":
                            parent_id = "NCBITaxon_" + row[7].split("NCBITAXON/")[1]
                            relationship = (node_id, parent_id)
                            edges.append(relationship)

                            # if len(node_parents) == 1: #
                            # Only consider concepts with 1 direct ancestor
                            child_to_parent[node_id] = parent_id

                        for synonym in synonyms:
                            synonym_to_id[synonym] = node_id

        # Create a MultiDiGraph object with only "is-a" relations
        # this will allow the further calculation of shorthest path lenght
        kb_graph = nx.DiGraph([edge for edge in edges])

        # Build id_to_info (KB-ID: (outdegree, indegree, num_descendants))
        for node in kb_graph.nodes:
            num_descendants = len(nx.descendants(kb_graph, node))

            id_to_info[node] = (
                kb_graph.out_degree(node),
                kb_graph.in_degree(node),
                num_descendants,
            )

        node_to_node = {}

        for edge in edges:
            if edge[0] in node_to_node:
                node_to_node[edge[0]].append(edge[1])

            else:
                node_to_node[edge[0]] = [edge[1]]

            if edge[1] in node_to_node:
                node_to_node[edge[1]].append(edge[0])

            else:
                node_to_node[edge[1]] = [edge[0]]

        self.name_to_id = name_to_id
        self.id_to_name = id_to_name
        self.id_to_info = id_to_info
        self.synonym_to_id = synonym_to_id
        self.edges = edges
        self.graph = kb_graph
        self.node_to_node = node_to_node

    def load_ncbi_gene(self):
        name_to_id = {}
        id_to_name = {}
        id_to_info = {}
        synonym_to_id = {}

        with open(f"{data_dir}ncbi_gene/All_Data.gene_info") as ncbi_gene:
            reader = csv.reader(ncbi_gene, delimiter="\t")
            row_count = int()

            for row in reader:
                row_count += 1

                if row_count > 7:
                    # Skip the header
                    gene_symbol = row[2]
                    gene_id = "NCBIGene_" + row[1]
                    synonyms = row[4].split("/")
                    description = row[8]
                    synonym_to_id[description] = gene_id

                    name_to_id[gene_symbol] = gene_id
                    id_to_name[gene_id] = gene_symbol

                    for synonym in synonyms:
                        if synonym != "-":
                            synonym_to_id[synonym] = gene_id

        edges = [("NCBIGene1", "NCBIGene2")]
        kb_graph = nx.DiGraph([edge for edge in edges])

        self.name_to_id = name_to_id
        self.id_to_name = id_to_name
        self.id_to_info = id_to_info
        self.synonym_to_id = synonym_to_id
        self.graph = kb_graph
        self.node_to_node = {}

    def load_txt(self):
        """Load knolwedge base from text files: terms.txt and edges.txt"""
        name_to_id = {}
        id_to_name = {}
        id_to_info = {}
        synonym_to_id = {}
        edges = []

        # import concept names
        with open(self.terms_filename, "r", encoding="utf-8") as in_file:
            data = in_file.readlines()

            for line in data:
                if line != "\n":
                    line_ = line.strip("\n").split("\t")
                    kb_id = line_[0]
                    name = line_[1]
                    name_to_id[name] = kb_id
                    id_to_name[kb_id] = name

                    if len(line_) == 3:
                        synonyms = line_[2].split(";")

                        for synonym in synonyms:
                            synonym_to_id[synonym] = kb_id

        # import relations between concepts
        edges = []

        with open(self.edges_filename, "r", encoding="utf-8") as in_file:
            data = in_file.readlines()

            for line in data:
                if line != "\n":
                    line_ = line.strip("\n").split("\t")
                    term1 = line_[0]
                    term2 = line_[1]
                    edges.append((term1, term2))

        kb_graph = nx.DiGraph([edge for edge in edges])

        # Build id_to_info (KB-ID: (outdegree, indegree, num_descendants))
        for node in kb_graph.nodes:
            num_descendants = len(nx.descendants(kb_graph, node))

            id_to_info[node] = (
                kb_graph.out_degree(node),
                kb_graph.in_degree(node),
                num_descendants,
            )

        node_to_node = {}

        for edge in edges:
            if edge[0] in node_to_node:
                node_to_node[edge[0]].append(edge[1])

            else:
                node_to_node[edge[0]] = [edge[1]]

            if edge[1] in node_to_node:
                node_to_node[edge[1]].append(edge[0])

            else:
                node_to_node[edge[1]] = [edge[0]]

        self.name_to_id = name_to_id
        self.id_to_name = id_to_name
        self.id_to_info = id_to_info
        self.synonym_to_id = synonym_to_id
        self.edges = edges
        self.node_to_node = node_to_node
        self.graph = kb_graph
