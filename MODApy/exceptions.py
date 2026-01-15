"""
Custom exceptions for MODApy.
"""

class MODApyError(Exception):
    """Base exception for all MODApy errors."""
    pass

class PatientNotFoundError(MODApyError):
    """Raised when a patient is not found in the database or directory."""
    pass

class PatientAlreadyExistsError(MODApyError):
    """Raised when attempting to add a patient that already exists in the database."""
    pass

class FileNotFoundError(MODApyError):
    """Raised when a required file is not found."""
    pass

class FileFormatError(MODApyError):
    """Raised when a file has an invalid format or extension."""
    pass

class DatabaseError(MODApyError):
    """Raised when database operations fail."""
    pass

class PanelError(MODApyError):
    """Raised when panel operations fail (e.g., no variants found for gene list)."""
    pass

class PipelineError(MODApyError):
    """Raised when pipeline execution fails."""
    pass

class VCFParseError(MODApyError):
    """Raised when VCF parsing fails."""
    pass
