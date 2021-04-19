
class StatsExtractionEnum:
    MOLECULAR_WEIGHTS = "mol_wts"
    NUMBER_OF_RINGS = "num_rings"
    NUMBER_OF_AROMATIC_RINGS = "num_aromatic_rings"
    NUMBER_OF_ATOMS = "num_atoms"
    HYDROGEN_BOND_DONORS = "hbond_donors"
    HYDROGEN_BOND_ACCEPTORS = "hbond_acceptors"
    HETERO_ATOM_RATIO = "hetero_atom_ratio"
    TOKEN_DISTRIBUTION = "token_distribution"
    TOKEN_ATOM_RATIO = "token_atom_ratio"
    NUMBER_OF_TOKENS = "num_tokens"
    NUMBER_OF_DECORATIONS = "num_decorations"

    REGEX_TOKENS = "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%\d\d|10|[0-9])"
    DECORATION_SEPARATOR_TOKEN = "|"
    ATTACHMENT_POINT_TOKEN = "*"

    SLICED_DATA_MODE = "sliced_data"
    SINGLE_COLUMN_DATA_MODE = "orig_data"

    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


