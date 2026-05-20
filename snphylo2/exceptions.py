"""
Custom exceptions for SNPhylo2.
"""


class SNPhylo2Error(Exception):
    """Base exception for all SNPhylo2 errors."""
    
    def __init__(self, message: str, details: dict | None = None):
        super().__init__(message)
        self.message = message
        self.details = details or {}


class ConfigurationError(SNPhylo2Error):
    """Raised when configuration is invalid or inconsistent."""
    pass


class InputError(SNPhylo2Error):
    """Raised when input files are invalid, missing, or malformed."""
    pass


class FormatError(SNPhylo2Error):
    """Raised when file format is not recognized or invalid."""
    pass


class FilterError(SNPhylo2Error):
    """Raised when filtering operations fail."""
    pass


class PruningError(SNPhylo2Error):
    """Raised when LD pruning fails."""
    pass


class TreeError(SNPhylo2Error):
    """Raised when tree building fails."""
    pass


class EngineError(SNPhylo2Error):
    """Raised when external tool execution fails."""
    
    def __init__(self, message: str, command: str, returncode: int, stderr: str = ""):
        super().__init__(message)
        self.command = command
        self.returncode = returncode
        self.stderr = stderr


class PopulationError(SNPhylo2Error):
    """Raised when population genetics analyses fail."""
    pass


class VisualizationError(SNPhylo2Error):
    """Raised when visualization generation fails."""
    pass


class ReportError(SNPhylo2Error):
    """Raised when report generation fails."""
    pass


class ValidationError(SNPhylo2Error):
    """Raised when validation checks fail."""
    pass
