"""Standalone LOFTEE-compatible transcript classification."""

from .classifier import classify
from .model import Classification, TranscriptContext

__all__ = ["Classification", "TranscriptContext", "classify"]
