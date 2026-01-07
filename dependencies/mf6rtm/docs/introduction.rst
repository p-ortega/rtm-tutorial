Introduction
============

Overview
--------

**mf6rtm** (MODFLOW 6 Reactive Transport Modeling) is a Python package that provides seamless integration between MODFLOW-6 and PHREEQC for simulating reactive transport processes in subsurface environments.

Reactive transport modeling is essential for understanding the intricate interplay between hydrogeological processes and chemical reactions in subsurface environments. mf6rtm bridges MODFLOW-6, the current version of the MODFLOW family of groundwater flow and transport codes, with PHREEQC, a versatile software for geochemical modeling.


What is Reactive Transport Modeling?
-------------------------------------

Reactive transport modeling combines:

* **Groundwater flow** - Movement of water through porous media
* **Solute transport** - Migration of dissolved components
* **Chemical reactions** - including:
  
  * Mineral dissolution and precipitation
  * Redox reactions
  * Ion exchange and sorption

Key Features
------------

Seamless Integration
~~~~~~~~~~~~~~~~~~~~

Through the integration facilitated by the MODFLOWAPI and PHREEQCRM APIs, mf6rtm provides a unified platform for modeling groundwater flow, solute transport, 
and chemical reactions within a single computational environment.

Uncertainty Analysis
~~~~~~~~~~~~~~~~~~~~

The code is designed to seamlessly integrate with PEST++ and PyEMU for:

* Uncertainty quantification
* Sensitivity analysis
* Parameter estimation
* Model calibration

This integration enables users to perform rigorous assessment of the impact of parameter and model uncertainties on reactive transport simulations.

Code Structure
------------

mf6rtm is organized into several subpackages, but there are two main components, that the user will interact with the most: solver and mup3d.

* **Solver**: This component manages the coupling between MODFLOW-6 and PHREEQC, handling data exchange, time-stepping, and overall simulation control. 
* **Mup3d**: This module acts as a pre- and post-processor for preparing input files for the reactive transport simulations, especially the chemistry inputs for PHREEQCRM (think FloPy for MODFLOW-6 + PHREEQC).

In addition to these, the user will require some knowledge and familiarity with Modflow 6 and FloPy to be able to set up and run mf6rtm simulations.

Why mf6rtm?
-----------

mf6rtm represents a significant advancement in reactive transport modeling by:

* **Integrating** state-of-the-art flow and geochemical codes
* **Providing** a Python-based, user-friendly interface
* **Enabling** uncertainty analysis and model calibration
* **Offering** flexibility for diverse hydrogeological scenarios

Getting Started
---------------

.. To get started with mf6rtm, see the :doc:`tutorials/index` (if available) 
or explore the :doc:`api/modules` documentation.

Installation
------------

Install mf6rtm using pip:

.. code-block:: bash

   pip install mf6rtm

For development installation:

.. code-block:: bash

   git clone https://github.com/p-ortega/mf6rtm.git
   cd mf6rtm
   pip install -e .
