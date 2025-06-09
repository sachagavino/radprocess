RAMSES AMR file formats
=======================

Default
-------

The default settings for the AMR data file formats is as follow::

    from pymses.sources.ramses.output import *
    RamsesOutput.amr_field_descrs_by_file = \
    {   "2D": {"hydro" : [ Scalar("rho", 0), Vector("vel", [1, 2]), Scalar("P", 3) ],
                "grav"  : [ Vector("g", [0, 1]) ]
            },
        "3D": {"hydro" : [ Scalar("rho", 0), Vector("vel", [1, 2, 3]), Scalar("P", 4) ],
             "grav"  : [ Vector("g", [0, 1, 2]) ]
            }
    }

which means that in the :file:`hydro_*.out*` files :
 * the first read variable corresponds to the scalar gas **density** field
 * the next 3 read variables corresponds to the gas 3D **velocity** field
 * the fifth read variable corresponds to the scalar gas **pressure** field

and in the :file:`grav_*.out*` files :
 * the 3 read variables corresponds to the 3D **gravitational acceleration** field


User-defined
------------

If you use a nD (:math:`n not equal to 3`) or a non-standard version of RAMSES, you might want to redefine this AMR file format to :
 * make additional tracers available to your reader
 * read nD (:math:`n not equal to 3`) data

::

    from pymses.sources.ramses.output import *
    # 2D data format
    RamsesOutput.amr_field_descrs_by_file = \
    {   "2D": {"hydro" : [ Scalar("rho", 0), Vector("vel", [1, 2]), Scalar("P", 3) ],
                "grav"  : [ Vector("g", [0, 1]) ]
            }
    }
    
    # Read additional tracers : metallicity, HCO abundancy
    RamsesOutput.amr_field_descrs_by_file = \
    {   "3D": {"hydro" : [ Scalar("rho", 0), Vector("vel", [1, 2, 3]), Scalar("P", 4), Scalar("Z", 5), Scalar("HCO", 6) ],\
             "grav"  : [ Vector("g", [0, 1, 2]) ]
            }
    }

To take into effect these settings, make sure you define them before any call to :meth:`~pymses.sources.ramses.output.RamsesOutput.amr_source`::

    from pymses.sources.ramses.output import *
    RamsesOutput.amr_field_descrs_by_file = \
    {   "3D": {"hydro" : [ Scalar("rho", 0), Vector("vel", [1, 2, 3]), Scalar("P", 4), Scalar("Z", 5)],\
             "grav"  : [ Vector("g", [0, 1, 2]) ]
            }
    }
    ro = RamsesOutput("/data/metal_simu/run001", 20)
    amr = ro.amr_source(["rho", "Z"])
    
User-defined output file
------------------------
In order to have user defined variables updated automatically by PyMSES when an output is opened, you can create a "pymses_field_descrs.py" in your current directory with this structure :

.. code-block:: python

    from pymses.sources.ramses import output
    self.amr_field_descrs_by_file = \
    	{"2D": {"hydro" : [ output.Scalar("rho", 0), output.Vector("vel", [1, 2, 3]),
    		output.Vector("Bl", [4,5,6]), output.Vector("Br", [7,8,9]),
    		output.Scalar("P", 10),output.Scalar("Z", 11)],
    			"grav"  : [ output.Vector("g", [0, 1, 2]) ]
    			},
    	"3D": {"hydro" : [ output.Scalar("rho", 0), output.Vector("vel", [1, 2, 3]),
    		output.Vector("Bl", [4,5,6]), output.Vector("Br", [7,8,9]),
    		output.Scalar("P", 10),output.Scalar("Z", 11)],
    			"grav"  : [ output.Vector("g", [0, 1, 2]) ]
    			}
    	}
    print "amr_field_descrs_by_file updated for MHD !"

