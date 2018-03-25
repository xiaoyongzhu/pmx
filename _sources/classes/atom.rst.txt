Atom
====

This class stores all properties (attributes) related to an atom.
When **pmx** reads a PDB file, each line starting with "ATOM" or "HETATM"
is converted to an :class:`Atom` instance. In addition, the atom class also has a number
of methods that can be used to calculate distances, angles and dihedrals.

API Reference
-------------

.. currentmodule:: pmx.atom

.. autoclass:: Atom

   .. rubric:: Methods

   .. autosummary::
      readPDBString
      dist
      dist2
