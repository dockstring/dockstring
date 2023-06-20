class DockstringWarning(Warning):
    """Base warning for Dockstring"""
    pass


class DockstringError(Exception):
    """Base exception for Dockstring"""
    pass


class CanonicalizationError(DockstringError):
    """Input SMILES cannot be canonicalized"""
    pass


class ParsingError(DockstringError):
    """Error raised when converting some input (file or string) to an RDKit Mol"""
    pass


class OutputError(DockstringError):
    """Error raised when writing an RDKit Mol to file or string"""
    pass


class SanityError(DockstringError):
    """Error raised when the input molecule doesn't fulfill certain criteria"""
    pass


class EmbeddingError(DockstringError):
    """Error raised during generation of 3D conformation"""
    pass


class StructureOptimizationError(DockstringError):
    """Error raised during structure optimization"""
    pass


class FormatConversionError(DockstringError):
    """Error raised when converting file formats"""
    pass


class ProtonationError(DockstringError):
    """Error raised when protonating a molecule"""
    pass


class PoseProcessingError(DockstringError):
    """Error raised during pose processing (e.g., assignment of stereochemistry)"""
    pass


class VinaError(DockstringError):
    """Error if AutoDock Vina fails"""
    pass


class DockingError(DockstringError):
    """Error raised if AutoDock Vina cannot find a pose"""
    pass
