r"""
This package is used for generating discontinuity factors after
receiving diffusion coefficients,cross-sections and interface currents from a transport code.

The purpose of discontinue factors is to compensate for error coming from the homogenization process.

Discontinuity factors are computed using compatible cross-sections, multiplication factor
and surface averaged currents. The computation process is as follows:

 - A transport code is used to compute the multiplication factor,
   the homogenized cross-sections of each cell and the surface averaged
   currents between different cells.
 - For each cell a homogenous diffusion problem is solved using the transport computed
   currents as boundary conditions.
 - For each cell :math:`c` and face :math:`f` the flux on the boundary :math:`\phi_{c,f}` is found.
 - The discontinuity factors are set to :math:`\frac{1}{\phi_{c,f}}`

"""

from dorban.rdfs.nem_rdf import compute_rdfs as compute_rdfs
