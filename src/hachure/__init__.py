"""Hachure engine: scalar/vector field -> Lehmann/Dufour-style strokes.

Pipeline: dem.fetch_copdem -> field.prepare -> trace.hachures -> export.*
The same engine is meant to run on NISAR-derived rasters (LOS displacement,
backscatter change) by passing them as the scalar field instead of elevation.
"""
