.. DORBAN documentation master file, created by
   sphinx-quickstart on Thu Sep 30 11:13:10 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
DORBAN
======

Welcome to DORBAN's documentation!

DORBAN is a python code used for modeling nuclear reactors and solving the various
forms of the diffusion approximation of the reactor transport equation.

DORBAN allows users to model reactors with varying geometries including Cartesian and Triangular
geometries, and allows solving different types of problems, including the direct and
adjoint flux problems as well as source problems, given both as a source inside the
core and as boundary conditions.

DOBRAN focuses on being readable, maintainable and modular rather on being fast.
DORBAN is written in python which allows high readability and maintainability, it also
allows the use of short and simple inputs to model complex system as shown in the
examples section.

==================================

.. toctree::
   :maxdepth: 2

   installation
   dorban
   examples
